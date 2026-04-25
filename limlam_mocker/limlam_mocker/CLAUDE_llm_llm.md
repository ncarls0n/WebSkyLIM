# `limlam_mocker/limlam_mocker/` — Architecture & Physics Reference

This is the inner Python package (imported as `import limlam_mocker as llm`). It implements the core simulation pipeline that converts a halo catalogue into a 3D line-intensity map and then into a power spectrum. The entry point in `limlam_mocker/lim_mocker.py` calls these modules in sequence.

---

## Pipeline Overview

```
halo catalogue (.npz / .h5)
        │
        ▼
load_halos.py          ← read file, compute RA/Dec/z, cull to map footprint
        │
        ▼
halos_to_luminosity.py ← assign L_line(M_halo, z) per halo
        │
        ▼
luminosity_to_map.py   ← bin (RA, Dec, ν_obs) with weights L into 3D cube
        │
        ▼
map_to_pspec.py        ← 3D FFT → spherically-averaged P(k)
```

---

## File-by-File Reference

### `__init__.py`
Wildcard-imports from every submodule so the public API is flat:
```python
from .tools               import *
from .load_halos          import *
from .halos_to_luminosity import *
from .luminosity_to_map   import *
from .map_to_pspec        import *
from .                    import debug
```
All functions documented below are accessible as `llm.<function>`.

---

### `debug.py`
One global flag: `verbose = False`. Set `llm.debug.verbose = True` before running to enable progress print statements throughout the pipeline.

---

### `tools.py`
Utility classes and functions shared across the pipeline.

**`empty_table`** — A bare Python class used as a struct to attach arbitrary attributes. Both `halos` and `map` instances throughout the code are `empty_table()` objects.

**`params_to_mapinst(params)`** — Converts a parameter object into a `map` instance (an `empty_table`) with all grid geometry pre-computed:
- Angular bins: `pix_binedges_x/y` in degrees, centered on the origin (RA=0, Dec=0).
- Frequency bins: `nu_binedges` from `nu_i` down to `nu_f` with step `dnu = (nu_i - nu_f)/nmaps`.
- Pixel solid angle: `Ompix = (pix_size_x × pix_size_y)` in steradians (used for brightness temperature conversion).
- Redshift bounds: derived as `z_i = nu_rest/nu_i - 1`, `z_f = nu_rest/nu_f - 1`.

**Cosmology helpers** (standalone, no astropy dependency):
- `hubble(z, h, omegam)` — H(z) = 100h √[Ω_M(1+z)³ + Ω_Λ] in km/s/Mpc.
- `drdz(z, h, omegam)` — c/H(z) in Mpc, the comoving line element dχ/dz.
- `redshift_to_chi(z, cosmo)` — Numerically integrates dχ/dz over z ∈ [0, 10] (20,000 steps) then interpolates to get χ(z).
- `chi_to_redshift(chi, cosmo)` — Inverse of above; interpolates z(χ).

**`timeme`** — Decorator that prints wall-clock seconds for a function call (unused in current pipeline, kept for profiling).

**`plot_results`** — Quick matplotlib visualization of the central-frequency slice and/or power spectrum.

---

### `load_halos.py`
Reads PeakPatch halo catalogues and applies a spatial/mass/redshift cut.

#### `load_peakpatch_catalogue(halo_info, filetype='.npz')`
Populates a `halos` struct with:

| Attribute | Description |
|---|---|
| `halos.M` | Halo mass in M_sun |
| `halos.x_pos`, `y_pos`, `z_pos` | Comoving Cartesian positions in Mpc |
| `halos.chi` | Comoving distance χ = √(x²+y²+z²) in Mpc |
| `halos.redshift` | Observed redshift (includes peculiar velocities for `.npz`) |
| `halos.zformation` | Formation redshift (`.npz` only) |
| `halos.ra`, `halos.dec` | Sky coordinates in degrees |
| `halos.vx`, `vy`, `vz` | Peculiar velocities in km/s (`.npz` only) |

**RA/Dec calculation**: `ra = atan2(-x, z) × 180/π`, `dec = arcsin(y/χ) × 180/π`. This projects the catalogue (which points along the +z axis) onto the sky.

**Two file formats:**
- `.npz` (default, PeakPatch standard): reads arrays directly from the archive; cosmology comes from `cosmo_header` dict stored in the file.
- `.h5` (Martine Lokken's catalogues): splits halos into centrals and satellites; redshift is computed numerically from χ using Planck15 cosmology since it isn't stored directly.

Sanity assertions: `max(M) < 1e17 M_sun`, `max(z) < 10`.

#### `load_peakpatch_catalogue_cosmo(halo_info)`
Extracts cosmological parameters (Ω_M, Ω_B, Ω_Λ, h, n_s, σ_8) from the `.npz` header into a `cosmo` struct. Asserts flatness (Ω_M + Ω_Λ = 1).

#### `cull_peakpatch_catalogue(halos, min_mass, max_mass, mapinst, haloType='all')`
Applies a combined boolean mask to all halo arrays:
```
M > min_mass
M < max_mass
|RA|  ≤ fov_x / 2
|Dec| ≤ fov_y / 2
z_i ≤ z ≤ z_f
```
`haloType` selects `'all'`, `'cen'` (centrals only), or `'sat'` (satellites only) — relevant only for `.h5` catalogues that separate the two populations.

---

### `halos_to_luminosity.py`
Assigns a line luminosity to each halo. Contains both legacy models and a modern dispatcher.

#### `Mhalo_to_Lline(halos, model_type, model_name, model_par, sigma_scatter, ...)` ← **main entry point used by `limlam.py`**
Dispatches to one of two modes:

**`model_type='ML'`** (mass-luminosity): Calls `getattr(mass_luminosity, model_name)(halos.M, model_par, halos.redshift)`. This routes to the `mass_luminosity.py` models in the repo root (`TonyLi`, `COMAP_Fid`, `Tolgay2026CO`, `SilvaCII`, etc.). Then optionally applies log-normal scatter. For `TonyLi`, an extra SFR-scatter term is applied with `pwr=1/alpha`.

**`model_type='LF'`** (luminosity function): Ignores halo mass entirely; draws N luminosities from the provided `dndL` PDF using inverse transform sampling (`rand_sample`). This matches the `luminosity_functions.py` workflow in the repo root.

#### Legacy models (kept but not called by default)

**`Mhalo_to_Lco_Li(halos, coeffs)`** — Tony Li et al. 2016 (arXiv:1503.08833). Chain:
1. M_halo → SFR via bivariate spline on the Behroozi+13a/b table (`sfr_behroozi_release.dat`). Axes: log₁₀(M/M_sun) × log₁₀(1+z).
2. SFR → L_IR = SFR × 10^10 / δ_MF (converts M_sun/yr → L_sun with a mass function correction).
3. L_IR → L_CO' = L_IR^(1/α) × 10^(-β/α). L_CO' is the CO line luminosity in observer units (K km/s pc²).
4. L_CO = 4.9×10⁻⁵ × L_CO' (conversion to L_sun).
5. Two rounds of log-normal scatter: one on SFR, one on L_CO.

Default coefficients: δ_MF=1, α=1.37, β=-1.74, σ_SFR=0.3 dex, σ_CO=0.3 dex.

**`Mhalo_to_Lco_Padmanabhan(halos, coeffs)`** — Padmanabhan 2017 (arXiv:1706.01471). Double-power-law with redshift-evolving parameters:
```
L'(M, z) = 2n(z) × M / [(M/M₁(z))^(-b(z)) + (M/M₁(z))^(y(z))]
L_CO = 4.9e-5 × L'
```
Each parameter (M₁, n, b, y) evolves as `p(z) = p₁₀ + p₁₁ × z/(1+z)`.

#### Scatter utilities

**`add_log_normal_scatter(data, dex, pwr=1.0, seed=None)`**
Draws multiplicative scalings from a log-normal with:
- σ_ln = dex × ln(10)  (converts dex → natural log units)
- μ_ln = −σ_ln²/2     (ensures `<scatter> = 1` in linear space, i.e., mean-preserving)

Returns `data × randscaling^pwr`. The `pwr` parameter is used in the TonyLi model for the SFR scatter, where pwr = 1/α propagates SFR uncertainty through the power-law L_IR relation.

**`add_variable_log_normal_scatter(data, z_halo, interp_table_file)`**
Reads σ(z) from a comma-separated table (columns: z, A, B, log C, log M*, σ) and applies per-halo scatter with redshift-interpolated σ. **Contains a hardcoded path `/home/njcarlson/...`** — must be given an explicit `interp_table_file` path outside the original HPC environment.

**`get_sfr_table(bad_extrapolation=False)`** — Loads `tables/sfr_behroozi_release.dat` (columns: 1+z, log M_halo, log SFR, log M_stellar) and builds a bilinear `RectBivariateSpline` on (log M, log(1+z)) axes. With `bad_extrapolation=True`, uses a degree-4 `SmoothBivariateSpline` to extend into the unphysical −1000 sentinel regions of the table.

**`rand_sample(x, PDF, N)`** — Inverse transform sampling: integrates the PDF to get the CDF, then interpolates CDF⁻¹(u) for N uniform random draws u ∈ [0,1].

---

### `luminosity_to_map.py`
Bins halo luminosities into the 3D intensity cube.

#### `Lco_to_Lmap(halos, map)`
1. Computes observed frequency for each halo: `ν_obs = ν_rest / (1 + z)`.
2. Calls `np.histogramdd` on `(RA, Dec, ν_obs)` with `weights=halos.Lco`.
3. The frequency axis is flipped before and after histogramming because `np.histogram` requires monotonically increasing bins but the frequency axis runs from high to low (high ν = low z = near end of survey).

**Output**: 3D array with shape `(npix_x, npix_y, nmaps)` in units of L_sun per voxel.

#### `Ngal_to_map(halos, map)`
Same as `Lco_to_Lmap` but with unit weights — produces a halo number count map. Used for cross-correlations.

#### `T_line(halos, map)`
Converts per-halo luminosity to brightness temperature in the Rayleigh-Jeans limit:

```
T = (c² / 2 k_B ν_obs²) × I_line
I_line = L / (4π D_L² dν Ω_pix)
D_L = χ(1+z)  [comoving distance × (1+z)]
```

The pre-computed conversion factor is 2.63083 μK per (L_sun Mpc⁻² GHz⁻¹ GHz² sr⁻¹). **Note**: this function exists but is not called in the current pipeline; `Lco_to_Lmap` stores luminosity directly and `limlam.py` does the L→T conversion separately using a CLT-based formula from `line_model.py`.

#### `save_maps(map)`
Saves the map cube and all grid metadata to a `.npz` file at `map.output_file`.

---

### `map_to_pspec.py`

#### `map_to_pspec(map, cosmo)`
Computes the spherically-averaged 3D power spectrum from the intensity map.

**Step 1 — Comoving coordinate conversion:**
- Frequency bin edges → redshift via `z = ν_rest/ν − 1`, then → comoving distance χ(z) using `redshift_to_chi`.
- Angular bin edges → comoving transverse distance: `x_co = θ_rad × χ_mean` (flat-sky approximation; assumes no curvature so transverse comoving distance = line-of-sight comoving distance at the field centre).

**Step 2 — 3D FFT:**
```
P_3D(k_vec) = |FFT(T) × dx × dy × dz|² / V_survey
```
Uses `np.fft.rfftn` (real FFT, exploiting the real-valued temperature cube). The normalization factor is the survey volume `V = Δx × Δy × Δz` in comoving Mpc³.

**Step 3 — Spherical averaging:**
- Constructs `|k| = √(k_x² + k_y² + k_z²)` grid.
- Sets bin width `dk = max(dk_x, dk_y, dk_z)`.
- `np.histogram` accumulates P_3D values and mode counts into uniform k bins from 0 to k_max.
- Returns `P(k) = Σ P_3D / N_modes`.

**Output dict**: `{k, Pk, nmodes}` where `k` and `Pk` are bin-centered arrays in units of Mpc⁻¹ and μK² Mpc³ respectively, and `nmodes` is the number of independent Fourier modes averaged in each bin (used for sample variance error estimates: `σ_P = P/√N_modes`).

---

### `mock_cmass_and_redmagic_gals.py`
**Standalone HPC script** — not part of the importable package. Hardcoded paths to CITA Niagara cluster data. Creates mock photometric galaxy catalogues (CMASS and redMaGiC) from PeakPatch halo catalogues for cross-correlation studies.

**CMASS mode**: Selects central halos with a Gaussian mass selection (log M_halo ~ N(12.78, 0.37)) matching the observed CMASS stellar mass distribution. Number of selected halos per redshift slice is set by scaling a small observed sky patch to full sky (`scale_up = 41253 deg² / box_area`).

**redMaGiC mode**: Selects both central and satellite halos in the mass range 2×10¹³/h < M < 3×10¹³/h with a uniform mass distribution per distance slice. Satellite fraction is set by `frac_cens = 0.5`. Number counts are again scaled to match the observed redMaGiC density.

**Dependencies**: `healpy`, `h5py`, `astropy`, a local `cosmology` module (not standard). Cannot be run outside the original HPC environment.

---

### `htl_v1.py`
Older copy of `halos_to_luminosity.py` kept for reference. Has the same `Mhalo_to_Lco_Li` and `Mhalo_to_Lco_Padmanabhan` functions but with the original (simpler) `get_sfr_table` that lacks the `bad_extrapolation` option. **Not imported by `__init__.py`** and not used in the current pipeline.

---

### `Ex2_hh.ipynb`
Jupyter notebook demonstrating a worked example using these modules. Kept in the package directory but not part of the importable API.

---

## Supporting Data Tables

Located at `limlam_mocker/tables/`:

| File | Used by | Description |
|---|---|---|
| `sfr_behroozi_release.dat` | `get_sfr_table()` | Behroozi+13a/b SFR(M_halo, z) on a grid of 1+z × log M; columns: 1+z, log M, log SFR, log M_stellar |
| `Tolgay2026_FIRE_CO.txt` | `mass_luminosity.Tolgay2026CO` (repo root) | Redshift-indexed double-power-law fit parameters from FIRE simulations; columns: z, A, B, log C, log M*, σ |

---

## Data Flow Summary

```
halos.M [M_sun]  +  halos.redshift
        │
        │  Mhalo_to_Lline (or legacy Mhalo_to_Lco_*)
        │  + log-normal scatter in dex
        ▼
halos.Lco [L_sun]  +  halos.ra, halos.dec
        │
        │  np.histogramdd with weights=Lco
        │  bins: (pix_binedges_x, pix_binedges_y, nu_binedges)
        ▼
map.maps [L_sun/voxel], shape (npix_x, npix_y, nmaps)
        │
        │  np.fft.rfftn + spherical averaging
        │  angular bins → comoving via χ_mean × θ_rad
        │  freq bins    → comoving via χ(nu_rest/nu - 1)
        ▼
k [Mpc⁻¹], Pk [μK² Mpc³], nmodes
```

## Notes for Navigation

- The `limlam.py` class in the repo root wraps this entire package. `LimLam.maps` calls `load_halos`, culls, calls `Mhalo_to_Lline`, then `Lco_to_Lmap`. `LimLam._Pk_dict` calls `map_to_pspec`.
- The L_sun/voxel → μK brightness temperature conversion in the simulation path is done in `limlam.py` using the `CLT` (conversion factor from line_model.py), not via `T_line` in this package.
- `add_variable_log_normal_scatter` has a **hardcoded absolute path** on line ~251 that must be overridden with an explicit `interp_table_file` argument.
- `mock_cmass_and_redmagic_gals.py` is a script, not a module — it executes at import time and will fail outside its original HPC environment.
