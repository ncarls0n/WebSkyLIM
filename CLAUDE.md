# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

WebSkyLIM is a Python toolkit for **Line Intensity Mapping (LIM)** — a cosmological technique that maps aggregate spectral-line emission (primarily CO and CII) from large cosmic volumes without resolving individual galaxies. It assigns line luminosities to dark matter halos, bins them into 3D intensity cubes, and compares to analytic forecasts for instruments like COMAP and CCAT-prime. Forked from Patrick Horlaville's CITA_LIM repo; builds on `limlam_mocker` (Stein, Chung, Breysse, et al.).

---

## Running and Testing

There is no build step or install. Run directly from the repo root with the virtual environment active:

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy scipy astropy hmf h5py multiprocess
```

Each core module contains **doctests** that run when invoked as `__main__`:

```bash
python lim.py
python line_model.py
python line_obs.py
python mass_luminosity.py
python luminosity_functions.py
```

There is no pytest, tox, or CI configuration.

Optional dependencies: `peakpatchtools` (PeakPatch-format halo catalogues, HPC-only), `camb` (non-linear matter power spectra via hmf).

---

## Architecture

The design is a **three-layer class hierarchy** with `cached_property` throughout for lazy, invalidation-aware computation.

```
lim()               # Factory in lim.py — primary entry point
  └── LineModel     # line_model.py — cosmology + astrophysics (z, H, CLT, dndM, LofM, Pk)
        └── LineObs # line_obs.py  — instrument/survey (Nch, Vvox, Pnoise, SNR, VID)
              └── LimLam  # limlam.py — simulation interface (halo catalogue → maps → Pk_sim)
```

**Caching system** (`_utils.py`): `cached_property` stores computed values in the instance dict and registers them in `_update_list`. Always use `m.update(**new_params)` to change parameters — direct attribute mutation bypasses cache invalidation.

**Simulation data flow in `LimLam`:**
1. `halo_info` — raw catalogue file read (`.npz`, `.h5`, or PeakPatch binary)
2. `halos` — culled to mass/redshift/sky bounds
3. `L_halos` — luminosities assigned via `mass_luminosity.py` ML models + log-normal scatter
4. `maps` — halos binned into 3D (RA, Dec, frequency) cube via `np.histogramdd`
5. `_Pk_dict` — FFT-based 3D power spectrum in comoving coordinates

---

## Module-by-Module Reference

### `lim.py` — Factory and Beam Utilities

**`lim(model_params, doObs, doSim, match_sim_cosmo)`** — The primary entry point. Reads `model_params` (either a string name of a dict in `params.py` or a dict directly), strips keys that are invalid for the chosen class tier using `remove_invalid_params`, then instantiates and returns the appropriate object:
- `doObs=False, doSim=False` → `LineModel`
- `doObs=True, doSim=False` → `LineObs`
- `doObs=True, doSim=True` → `LimLam`

If `match_sim_cosmo=True` and the simulation cosmology differs from the analytic model, `set_sim_cosmo` is called to align them.

**`pix_res(m, beam_width, map_dim_pix, map_dim_deg)`** — Converts a beam FWHM to a pixel radius for `scipy.ndimage.gaussian_filter`. Returns `round(beam_width_arcmin / pixel_size_arcmin)`. Defaults to `m.beam_width`, `m.Nside`, `m.Omega_field` if not specified.

**`gaussian_beam(cube, m, ...)`** — Convenience wrapper: calls `pix_res` then applies `gaussian_filter` to the 3D intensity cube.

---

### `line_model.py` — Core Analytic Astrophysics (`LineModel`)

The base class. All properties are `@cached_property`. Uses `hmf.MassFunction` (Steven Murray's package) for halo abundances and `astropy.cosmology` for distances.

**Cosmology and target redshift:**
- `z` — Emission redshift: `z = nu/nuObs − 1`
- `H` — Hubble parameter H(z) from astropy cosmology
- `CLT` — Conversion factor from luminosity density to brightness temperature (or Jy/sr if `do_Jysr=True`):
  - Temperature mode: `CLT = c³(1+z)² / (8π k_B ν³ H)` in uK Mpc³/L_sun
  - Jy/sr mode: `CLT = c / (4π ν H sr)` in Jy Mpc³/(L_sun sr)
- `h` — `hmf.MassFunction` instance (not cleared by `update()` — cosmology changes go through hmf's own updater)

**Mass/luminosity grids:**
- `M` — Logarithmically spaced mass grid from `Mmin` to `Mmax` with `nM` points
- `L` — Logarithmically spaced luminosity grid from `Lmin` to `Lmax` with `nL` points
- `k`, `k_edge`, `dk` — Wavenumber grid (log or linear depending on `linear_k`)
- `dndM` — Halo mass function in Mpc⁻³ Msun⁻¹, interpolated from hmf output (converts from 1/h units)
- `sigmaM` — Mass variance σ(M) from hmf, used for bias calculation
- `LofM` — L(M) evaluated at the mass grid: calls `getattr(ml, model_name)(M, model_par, z)` for ML models, or a linear-in-M proxy for LF models (used only for bias weighting)
- `dndL` — Luminosity function in Mpc⁻³ L_sun⁻¹. For LF models: calls `getattr(lf, model_name)(L, model_par)`. For ML models: inverts L(M) numerically assuming monotonicity and applies scatter convolution.

**Halo bias:**
- `bofM` — Tinker et al. 2010 bias formula: `b(M) = 1 + (ν² − 1)/δ_c` where `ν = δ_c/σ(M)` and `δ_c` is the collapse threshold from hmf
- `bavg` — Luminosity-weighted average bias: `∫ L(M) b(M) dn/dM dM / ∫ L(M) dn/dM dM`

**NFW profile (one-halo term):**
- `ft_NFW` — Fourier transform of the NFW density profile on the (k, M) grid, using the Wechsler et al. 2002 concentration–mass relation `c = 4.1 / (a_c(1+z))`. Returns the normalized profile `u(k|M)`.

**Signal statistics:**
- `nbar` — Mean galaxy number density; from `∫ dn/dL dL` (LF) or `∫ dn/dM dM` (ML)
- `Tmean` — Mean brightness temperature: `CLT × ∫ L dn/dL dL` or `CLT × ∫ L(M) dn/dM dM`. For TonyLi, a scatter correction factor `exp((α⁻² − α⁻¹) σ²_SFR ln²10 / 2)` is applied to account for the non-mean-preserving SFR scatter.
- `Pshot` — Shot noise power: `CLT² × ∫ L² dn/dL dL`. For ML models, scatter inflates this by `exp(σ²_scatter ln²10)`.
- `Pk_twohalo` — Two-halo clustering power: `(T_mean × b_avg)² × P_m(k)`. If `do_onehalo=True`, uses the NFW-weighted integral `[CLT × ∫ L(M) b(M) u(k|M) dn/dM dM]²`.
- `Pk_onehalo` — One-halo term: `CLT² × ∫ L²(M) u²(k|M) dn/dM dM`. Zero if `do_onehalo=False`.
- `Pk` = `Pk_twohalo + Pk_onehalo + Pshot`
- `Pm` — Matter power spectrum from hmf (linear, or non-linear if `nonlinear_Pm=True`)

**1D power spectra:** `P1D_clust`, `P1D_shot`, `P1D` — line-of-sight (radial) power spectrum, integrated as `P_{1D} = ∫_{k>k_perp} k P(k) dk / (2π)`.

**`update(**new_params)`** — The only safe way to change parameters. Clears all cached properties in `_update_list`, sets new attribute values, and updates the hmf object for cosmology/redshift changes.

**`reset()`** — Calls `update` with the original `__init__` inputs.

---

### `line_obs.py` — Instrument Layer (`LineObs(LineModel)`)

Inherits all of `LineModel`. Adds survey geometry, instrument noise, and the Voxel Intensity Distribution (VID). Default parameters match COMAP1.

**Survey geometry:**
- `Nch` — Number of frequency channels: `round(Delta_nu / dnu)`
- `Nside` — Pixels per map side: `round(sqrt(Omega_field) / beam_width)`. Pixel size = one beam FWHM.
- `Nvox = Npix × Nch` — Total voxels
- `Vfield` — Comoving survey volume: `r₀² × Omega_field × c × Delta_nu × (1+z)² / (H × nu)` in Mpc³
- `Vvox = Vfield / Nvox` — Comoving voxel volume

**Instrument noise:**
- `sigma_N` — RMS noise per voxel:
  - Temperature mode: `T_sys / sqrt(N_feeds × dnu × t_pix)` in μK
  - Jy/sr mode: `NEFD / (Omega_beam × sqrt(N_feeds × t_pix))` where `Omega_beam = π FWHM² / (4 ln2)`
- `Pnoise = sigma_N² × Vvox` — Noise power spectrum amplitude
- `sigma_par`, `sigma_perp` — Physical resolution limits along and across the line of sight, in Mpc. `sigma_par = c × dnu_total × (1+z) / (H × nuObs)` where `dnu_total` includes intrinsic line width.
- `Wk` — Beam/channel window function in Fourier space: `exp(−k² σ_perp²) × ∫₀¹ exp(−k² (σ_par² − σ_perp²) μ²) dμ`

**Power spectrum errors:**
- `Nmodes = k² dk V_field N_field / (4π²)` — Number of independent Fourier modes per k bin
- `sk_CV = Pk / sqrt(Nmodes)` — Sample variance error
- `sk_N = Pnoise / sqrt(Nmodes)` — Noise error
- `sk = sk_CV + sk_N` — Total error (conservative sum, not in quadrature)
- `SNR` — Total signal-to-noise: `sqrt(sum over k of [Pk_obs / sk_obs]²)`

**Angular power spectrum (`Cl`):**
- `_Cl_exact` — Full integral: `(2/π) ∫ k² P(k) [∫ j_l(kr) dr]² dk` using parallelized (`multiprocess.Pool`) evaluation of spherical Bessel functions at each multipole l
- `_Cl_limber` — Limber approximation: `∫ H(z) P(l/r) / [c r² Δz²] dz`, cheaper but less accurate at low l
- `Cl` — Dispatches to exact or Limber based on `do_limber`
- `Clnoise` — Angular noise power; `sl`, `sl_CV`, `sl_N` — errors per multipole; `SNR_2D` — angular SNR

**Voxel Intensity Distribution (VID):**
The VID is the predicted histogram of voxel intensities, used as a non-Gaussian statistical estimator beyond the power spectrum.

- `Nbar = nbar × Vvox` — Mean number of emitters per voxel
- `PofN` — Probability of N galaxies in a voxel: convolves a log-normal density PDF (parametrized by `sigma_G`) with Poisson statistics: `P(N|μ̄) = ∫ Poisson(N|μ) × LogNormal(μ|μ̄, σ_G) dμ`
- `P1` — Single-emitter intensity PDF: `P1(T) = dn/dL(T / X_LT) / (n̄ × X_LT)` where `X_LT = CLT / Vvox`
- `PT` — Full voxel intensity PDF: convolution of P1 with itself N times, weighted by PofN. Fast mode (`do_fast_VID=True`): `PT = FFT⁻¹(Σ_N PofN[N] × FFT(P1)^N)`. Slow mode: direct convolution.
- `PT_zero = PofN[0]` — Probability of empty voxels (delta function at T=0, tracked separately)
- `Bi` — Predicted histogram bin counts: `Bi = Nvox × ∫_{bin i} PT dT`

---

### `limlam.py` — Simulation Interface (`LimLam(LineObs)`)

Wraps `limlam_mocker` to generate simulated maps and compares them against analytic predictions.

**Catalogue loading (`halo_info`):**
- `.h5` files: opened with `h5py.File` (Martine Lokken's centrals+satellites catalogues)
- `.npz` files: loaded with `np.load(..., allow_pickle=True)` (standard PeakPatch output)
- Other paths: loaded with the `peakpatchtools` package (HPC-only, reads PeakPatch binary `.pksc` format and converts to the standard dict structure)

**`limlam_params`** — Translates `LineObs` parameters into the `empty_table` struct expected by `limlam_mocker.params_to_mapinst`: `Nch→nmaps`, `Omega_field→fov_x/fov_y`, `Nside→npix_x/npix_y`, `nuObs±Delta_nu/2 → nu_i/nu_f`.

**`halos`** — Calls `llm.load_peakpatch_catalogue` then `llm.cull_peakpatch_catalogue` to cut to `[Mmin, Mmax]` and the map footprint. `haloType` controls whether to use all/centrals/satellites for `.h5` catalogues.

**`L_halos`** — Calls `llm.Mhalo_to_Lline` with `model_type`, `model_name`, `model_par`, `sigma_scatter`, and `scatter_seed`. For LF models, passes the analytic `L` and `dndL` arrays for inverse-transform luminosity sampling.

**`maps`** — Calls `llm.Lco_to_Lmap` (which runs `np.histogramdd`), then converts from L_sun/voxel to brightness temperature by multiplying by `X_LT = CLT / Vvox`.

**`noise_added_map`** — Applies a 1-pixel Gaussian beam smoothing across RA/Dec (not frequency), then adds Gaussian white noise with `sigma = sigma_N`.

**`_Pk_dict`** — Calls `llm.map_to_pspec`, which FFTs the map, converts angular bins to comoving coordinates via χ(ν_rest/ν − 1), and spherically averages. Returns `{k, Pk, nmodes}`.

**Simulation diagnostics:**
- `Tmean_sim = maps.mean()` — Mean map intensity
- `Tmean_halos = CLT × sum(L_halos) / Vfield` — Mean intensity from direct halo sum
- `Pshot_halos = CLT² × sum(L_halos²) / Vfield` — Shot noise from direct halo sum
- `dndL_sim`, `dndM_sim` — Simulated luminosity/mass functions (histogrammed)
- `Bi_sim` — Simulated VID histogram (from `np.histogram` on `maps`)
- `PofN_sim` — Simulated galaxy number count PDF (from `Ngal_maps`)

**`set_sim_cosmo(limlam)`** (module-level function) — Reads the cosmology from the halo catalogue and updates the `lim` object's `cosmo_model` to match, constructing an `astropy.FlatLambdaCDM` from the catalogue's `h`, `Omega_M`, `Omega_B`.

---

### `mass_luminosity.py` — L(M) Model Library

All functions have the signature `f(Mvec, MLpar, z)` where `Mvec` is in M_sun with astropy units and return L in L_sun with units. Called via `getattr(ml, model_name)(...)` in `LineModel`.

| Function | Line | Physics |
|---|---|---|
| `MassPow` | any | Simple power law: `L = A × (M/M_sun)^b` |
| `DblPwr` | any | Double power law with redshift evolution: `L = A × 10^(b1 z) × (M/10⁸)^b2 × (1 + M/M*)^b3` |
| `COMAP_Fid` | CO(1-0) | COMAP fiducial: `L = C / [(M/Ms)^A + (M/Ms)^B]` (arXiv:2111.05931, Table 1 "UM+COLDz+COPSS") |
| `Tolgay2026CO` | CO(1-0) | FIRE-calibrated double power law; reads `limlam_mocker/tables/Tolgay2026_FIRE_CO.txt`; interpolates A, B, log C, log M*, σ vs z; converts L_CO' (K km/s pc²) to L_sun via `L = L' × (8π k_B ν³ / c³)` with ν = 115.2712 GHz |
| `TonyLi` | CO(1-0) | Li et al. 2016 (arXiv:1503.08833): M→SFR (Behroozi+13 table)→L_IR→L_CO' via `L_CO' = (L_IR/10^β)^(1/α)`; `L_CO = 4.9e-5 L_CO'`. Two scatter terms: σ_SFR on SFR (non-mean-preserving) and σ_LCO on L_CO |
| `SilvaCII` | [CII] | Silva et al. 2015: `log L_CII = a × log(SFR) + b`. SFR from Behroozi+13 (via `Behroozi_SFR`). Uses `sfr_release.dat`. |
| `LichenCII` | [CII] | HI-mass based: `M_HI(M) = M0 × (M/Mmin)^α_HI × exp(−(Mmin/M)^0.35)`; `L_CII = α_LCII × M_HI` |
| `LichenCII_v3` | [CII] | Extension of LichenCII with additional bias parameters `alpha0`, `gamma0` and log-normal scatter `zdex` on the luminosity |
| `SilvaCO` | CO | Silva et al. CO model; piecewise power law in M with 5 segments |
| `Padmanabhan_CII` | [CII] | Padmanabhan 2018: `L = F(z) × (M/M1)^β × exp(−N1/M)` with redshift-scaling function F(z) |
| `Be13_Pwr` | any | Power law between L and SFR: `L = K × SFR(M,z)^γ`; SFR from Behroozi+13 |
| `CII_Metallicity` | [CII] | Collisional excitation model from Gong 2012 / Silva 2015; depends on electron temperature and density |
| `MassPow` | any | Bias-proxy model (linear in M); used internally by `LineModel` for LF-mode bias weighting |
| `Behroozi_SFR` | — | Helper: interpolates `sfr_release.dat` bivariate spline at given (M, z) |
| `Silva_SFR` | — | Helper: evaluates Silva+15 double-power-law SFR(M,z) from `Silva15_SFR_params.dat` (dormant) |

**`Tolgay2026CO` detail**: Uses `pathlib.Path(__file__).resolve().parent` to locate the table, so it is path-independent. Scatter is applied by drawing directly from `np.random.normal(lnL_CO', σ_ln)` — i.e., L_CO' is the mean of the log-normal, not L_CO'. The conversion constant `C = 8π k_B ν³_CO / c³ = 4.906e-5 L_sun/(K km/s pc²)` is derived from fundamental constants.

---

### `luminosity_functions.py` — dn/dL Parameterizations

All functions have signature `f(Lvec, LFparams)` and return dn/dL in Mpc⁻³ L_sun⁻¹. Called via `getattr(lf, model_name)(...)` in `LineModel`.

| Function | Form |
|---|---|
| `Sch` | Schechter: `φ* (L/L*)^α exp(−L/L*)` with hard low-L cutoff at Lmin |
| `SchCut` | Schechter with exponential double cutoff: `φ* (L/L*)^α exp(−L/L* − Lmin/L)` |
| `FromFile` | Reads two-column file (L in L_sun, dn/dL in Mpc⁻³ L_sun⁻¹); interpolates with `interp1d` |

---

### `params.py` — Parameter Preset Dicts

Flat Python dicts consumed by `lim()`. All astropy units are embedded. Key presets:

| Name | Line | Instrument | Notes |
|---|---|---|---|
| `default_par` | CO(1-0) | COMAP-like | Default LF model (SchCut from Breysse+2017) |
| `TonyLi_PhI` | CO(1-0) | COMAP Phase I (Ku band, 15 GHz) | z ≈ 6.7, 4 fields |
| `TonyLi_PhII` | CO(1-0) | COMAP Phase II (Ka band, 26 GHz) | z ≈ 3.4, 538 deg² |
| `COMAP_Fid` | CO(1-0) | COMAP Phase II | Fiducial double-power-law L(M) from arXiv:2111.05931 |
| `Silva_m1/m2/m3/m4_z6_CCATp` | [CII] | CCAT-prime (250 GHz, z ≈ 6.6) | Four Silva+15 model variants; NEFD-based noise |
| `Lichen`, `Lichen_v2/v3/v4` | [CII] | CCAT-prime (270 GHz, z ≈ 6) | LichenCII model evolution; v4 is the production version used in Horlaville+2023 |
| `Hamsa_CII` | [CII] | COMAP-like | Padmanabhan+2018 CII model |
| `Silva_CO32_StageII` | CO(2-1) or CO(3-2) | CMB-S4-era | Far-future sensitivity assumptions |

Many presets have hardcoded `catalogue_file` paths to CITA HPC systems. Override locally before use:
```python
m = lim('COMAP_Fid', catalogue_file='/local/path/to/catalogue.npz', doSim=True)
```

---

### `_utils.py` — Shared Utilities

**`cached_property`** — Django-style descriptor. On first access, calls the method, stores the result in `instance.__dict__[name]` (bypassing the descriptor for all future accesses), and appends the property name to `instance._update_list`. `update()` uses that list to `delattr` cached values.

**`check_params(input_params, default_params)`** — Type and unit validation at `__init__` time. Raises `AttributeError` for unknown keys, `TypeError` for wrong types or incompatible units.

**`check_model(model_type, model_name)`** — Verifies that `model_name` exists as a function in the correct module (`mass_luminosity` for `'ML'`, `luminosity_functions` for `'LF'`).

**Array utilities:** `ulogspace`, `ulinspace` — astropy-unit-aware versions of `np.logspace`/`np.linspace`. `uinterp1d` — unit-aware wrapper around `scipy.interpolate.interp1d`.

**Matrix utilities:** `rebin_map`, `bin_ndarray` — spatial rebinning of 3D maps by averaging or summing. `cov_to_corr`, `plotCov`, `conc22` — covariance matrix tools.

**`LtoLprime(L, nu)`** — Converts L_sun to observer-frame CO luminosity units K km/s pc²: `L'= L / (4.9e-5) × (ν/115 GHz)^-3`.

**`WF` class** — Wiener filter: given data vector `d`, noise covariance `N`, and signal covariance `S`, computes filter `F = S(S+N)⁻¹`, filtered signal `s = Fd`, and signal covariance `C_s = FN`.

**`MidpointNormalize`** — Matplotlib colormap normalizer that maps a specified midpoint to 0.5 (useful for diverging colormaps centered on zero).

---

### `_vid_tools.py`

Helper functions for `LineObs` VID calculations. Contains `lognormal_Pmu`, `binedge_to_binctr`, `binctr_to_binedge_log`, `pdf_to_histogram`, and `conv_parallel` (used by the slow VID path).

---

### `stacking_utils.py` — Halo-Stacking Core Library

Functions for extracting and averaging intensity map cutouts centred on known halo positions. Used by all `stacking_*.py` scripts.

**`halo_centpix(lim_obj, halo_xpos, halo_ypos, halo_zpos)`** — Converts halo sky coordinates (RA, Dec, redshift) to pixel indices on the LIM map using `np.argmin(|halo_pos − map_grid|)` for each axis.

**`halo_map(lim_obj, n, ...)`** — Generates an `n×n` pixel cutout index grid centred on each halo's pixel position. Marks out-of-bounds pixels as `NaN`.

**`inbound_halos(...)`** — Returns the subset of halos whose entire `n×n` cutout lies within the map boundaries.

**`lum_hod(lim_obj, n, ...)`** — Main stacking function: extracts an `n×n` cutout from `lim_obj.maps` (pure signal) and `lim_obj.noise_added_map` at each halo position. Returns arrays of shape `(n_halos, n, n)`.

**`lum_hod_corr(...)`** — Variant that applies a signal scaling factor (`sig_scale = 1.23`) and re-generates noise internally via `own_noise`, rather than using `noise_added_map`.

**`own_noise(m, sign, sigm, pix=1.)`** — Applies a per-slice Gaussian beam (`gaussian_filter([pix, pix, 0])`), then adds Gaussian white noise with RMS = `sigm` converted to Jy/sr.

**`n_mh(m_h, logm_min, sigma_logm, alpha, dc)`** — HOD (Halo Occupation Distribution): returns expected number of galaxies per halo. Centrals: `N_c = (1/2)[1 + erf((log M − log M_min)/σ)]`. Satellites: `N_s = N_c × [(M − M_0)/M_1']^α`. Total: `N = dc × (N_c + N_s)` where `dc` is the duty cycle. `M_1' = 10^(1.18 log M_min − 1.28)`, `M_0 = 10^(0.76 log M_1' + 2.3)` (Harikane et al. 2016 LBG HOD fit).

**`pix_res(beam_width, map_dim_deg, map_dim_pix)`** — Local version of the `lim.py` function; same calculation but takes arguments directly instead of from a `lim` object.

---

### `stacking_params.py` — Stacking Configuration

Shared setup imported by all `stacking_*.py` scripts. Instantiates the `lim_sim` object with `Lichen_v4` parameters (CCAT-prime CII at nuObs=270 GHz), defines:
- `t_obs = 2000 hr` — total observation time
- `mass_cut = 2e10 M_sun` — minimum halo mass for stacking
- `err = 0.1` — redshift tolerance window for halo selection (|z_halo − z_map| < err)
- `n = 50` — cutout size in pixels
- `z_sel = 6` — target redshift slice
- `beam_width = 50 arcsec` — CCAT-prime 50'' beam for [CII]
- `log_m_min = 11.03`, `sigma_log_m = 0.2`, `alph = 1`, `dutyc = 0.6` — HOD parameters for LBGs at z=5.9 (Harikane+2016)

---

### `stacking_*.py` — Stacking Analysis Scripts

All scripts import from `stacking_params.py` (which instantiates the `lim_sim` object) and run HPC-dependent analysis. All have hardcoded catalogue paths.

| Script | What it stacks | Method |
|---|---|---|
| `stacking.py` | Halos above `mass_cut` in a single lightcone | `lum_hod`; plots pure and noisy stacked maps |
| `stacking_alt.py` | Same, using alternative `lum` function call | `lum` (older API); single lightcone |
| `stacking_many.py` | Ensemble average over all lightcones in a directory | Loops over `lc_paths`; saves each noisy stack as `.npy` |
| `stacking_many_corr.py` | Same but with signal correction (`lum_hod_corr`) | `sig_scale = 1.23` applied before noise |
| `stacking_LBGs.py` | HOD-selected LBG-like galaxies from halos | Draws Poisson realisations from `n_mh`; weighted average |
| `stacking_many_LBGs.py` | Ensemble HOD-weighted stacking over all lightcones | Combines `stacking_many` + HOD weighting |
| `stacking_many_LBGs_corr.py` | Same with signal correction | `lum_hod_corr` + HOD weighting |

---

### Root Data Files

| File | Used by | Description |
|---|---|---|
| `sfr_release.dat` | `mass_luminosity.TonyLi`, `SilvaCII`, `COMAP_Fid`, `Be13_Pwr` | Behroozi+13 SFR(M, z) grid; same data as `limlam_mocker/tables/sfr_behroozi_release.dat` but without the 2-line header |
| `sfr_reinterp.dat` | `LichenCII_v2/v3` (via `BehrooziFile` key) | Re-interpolated Behroozi table on a finer grid; referenced as `'BehrooziFile':'sfr_reinterp.dat'` in Lichen presets |

---

## Parameter System

`params.py` contains preset dicts (e.g., `default_par`, `COMAP_Fid`, `Tolgay2026CO`, `Silva_m1-m4_z6_CCATp`, `Lichen`). Many presets hard-code `catalogue_file` paths for HPC systems — these must be overridden locally:

```python
m = lim(model_params='COMAP_Fid', catalogue_file='/local/path/to/catalogue.npz', doSim=True)
```

---

## Data Directories

| Directory | Reference doc | Contents |
|---|---|---|
| `fire_galaxies/` | [`fire_galaxies/CLAUDE_fire_galaxies.md`](fire_galaxies/CLAUDE_fire_galaxies.md) | FIRE simulation CII+IR luminosities (z=6, z=8) for model calibration |
| `sfr_tables/` | [`sfr_tables/CLAUDE_sfr_tables.md`](sfr_tables/CLAUDE_sfr_tables.md) | `sfr_release.dat` (Behroozi+13 SFR grid), `Silva15_SFR_params.dat` (Silva+15 SFR double power law) |
| `readmes/` | [`readmes/CLAUDE_readmes.md`](readmes/CLAUDE_readmes.md) | Tutorial notebooks by Horlaville: `limlam_mocker` intro and relative entropy (PRE) analysis |
| `limlam_mocker/` | [`limlam_mocker/CLAUDE_llm_pkg.md`](limlam_mocker/CLAUDE_llm_pkg.md) | Vendored simulation package root — `lim_mocker.py`, `params.py`, subdirectory index |
| `limlam_mocker/limlam_mocker/` | [`limlam_mocker/limlam_mocker/CLAUDE_llm_llm.md`](limlam_mocker/limlam_mocker/CLAUDE_llm_llm.md) | Inner package: halo loading, luminosity assignment, map binning, power spectrum |
| `limlam_mocker/catalogues/` | [`limlam_mocker/catalogues/CLAUDE_catalogues.md`](limlam_mocker/catalogues/CLAUDE_catalogues.md) | `.npz` catalogue format, `split_halo_catalogue.py`, COMAP mid-z metadata |
| `limlam_mocker/tables/` | [`limlam_mocker/tables/CLAUDE_tables.md`](limlam_mocker/tables/CLAUDE_tables.md) | `sfr_behroozi_release.dat`, `Tolgay2026_FIRE_CO.txt` |

---

## Active Development Notes

- **`Tolgay2026CO`** in `mass_luminosity.py` (~line 112) is the newest model under active development. It reads `limlam_mocker/tables/Tolgay2026_FIRE_CO.txt` and uses redshift-interpolated double-power-law fits from FIRE simulations. The inner helper functions `logL_CO_prime_z_lt_2` / `z_gt_1` / `interp` use a separate (older) code path that only activates when `'interp_table_file'` is absent from `MLpar`; the table-based path (default) is the active one.
- `add_variable_log_normal_scatter` in `limlam_mocker/limlam_mocker/halos_to_luminosity.py` contains a hardcoded absolute path (`/home/njcarlson/...`) that will fail outside the original HPC environment.
- `main.py` is almost entirely commented out; `runner.py` requires `peakpatchtools` (HPC-only).
- The three Behroozi SFR table copies (`sfr_release.dat`, `sfr_reinterp.dat` at repo root; `limlam_mocker/tables/sfr_behroozi_release.dat`) must be kept consistent if the table is ever updated.
- `sfr_release.dat` and `Silva15_SFR_params.dat` are loaded by bare filename (no path), so `mass_luminosity.py` must be called from the repo root working directory.
