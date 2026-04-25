# `limlam_mocker/` — Package Root Reference

This directory is the top-level of the vendored `limlam_mocker` package — the original George Stein / Dongwoo Chung simulation library, adapted for use as a sub-package of WebSkyLIM. It is not installed; it is imported directly because the repo root is on the Python path.

```
limlam_mocker/            ← this directory (package root)
├── lim_mocker.py         ← standalone entry-point script (original usage)
├── params.py             ← parameter file for lim_mocker.py
├── limlam_mocker/        ← inner Python package (importable as `limlam_mocker`)
│   ├── CLAUDE_llm_llm.md ← architecture reference for the inner package
│   ├── __init__.py
│   ├── tools.py
│   ├── load_halos.py
│   ├── halos_to_luminosity.py
│   ├── luminosity_to_map.py
│   └── map_to_pspec.py
├── catalogues/           ← halo catalogue data and pre-processing tools
│   ├── CLAUDE_catalogues.md
│   ├── split_halo_catalogue.py
│   └── COMAP_z4.77-6.21_700Mpc/simdata_info.txt
└── tables/               ← tabulated astrophysical data
    ├── CLAUDE_tables.md
    ├── sfr_behroozi_release.dat
    └── Tolgay2026_FIRE_CO.txt
```

---

## `lim_mocker.py` — Original Standalone Entry Point

The original top-level script showing the canonical usage of `limlam_mocker` as a standalone tool. **Not used by the WebSkyLIM class hierarchy** (`lim.py` / `LimLam` in the repo root take this role), but preserved as a reference implementation and minimal working example.

```python
import limlam_mocker as llm
import params as params

llm.debug.verbose = True
mapinst   = llm.params_to_mapinst(params)           # build map grid from params
halos, cosmo = llm.load_peakpatch_catalogue(params.halo_catalogue_file)
halos     = llm.cull_peakpatch_catalogue(halos, params.min_mass, mapinst)
halos.Lco = llm.Mhalo_to_Lco(halos, params.model, params.coeffs)  # legacy API
mapinst.maps = llm.Lco_to_map(halos, mapinst)
llm.save_maps(mapinst)
k, Pk, Nmodes = llm.map_to_pspec(mapinst, cosmo)
Pk_sampleerr  = Pk / np.sqrt(Nmodes)
llm.plot_results(mapinst, k, Pk, Pk_sampleerr, params)
```

Note that this calls `Mhalo_to_Lco` — the legacy dispatcher that routes to `Mhalo_to_Lco_Li` or `Mhalo_to_Lco_Padmanabhan`. The WebSkyLIM `LimLam` class calls `Mhalo_to_Lline` instead, which routes through the `mass_luminosity.py` model library.

---

## `params.py` — Default Parameter File for `lim_mocker.py`

Flat Python module (no classes); its attributes are read directly by `lim_mocker.py` and `params_to_mapinst`. **This `params.py` is separate from the `params.py` in the repo root**, which holds WebSkyLIM parameter presets.

| Parameter | Default | Description |
|---|---|---|
| `model` | `'Li'` | L(M) model name (`'Li'` or `'Padmanabhan'`) |
| `coeffs` | `None` | Model coefficients; `None` uses the published defaults |
| `halo_catalogue_file` | `'catalogues/peakpatchcatalogue_1pt4deg_z2pt4-3pt4.npz'` | Path to `.npz` catalogue |
| `min_mass` | `2.5e10` | Minimum halo mass in M_sun |
| `nu_rest` | `115.27` | CO(1-0) rest frequency in GHz |
| `nu_i` | `34.0` | High-frequency (low-z) edge of survey band in GHz |
| `nu_f` | `26.0` | Low-frequency (high-z) edge of survey band in GHz |
| `nmaps` | `100` | Number of frequency channels |
| `fov_x`, `fov_y` | `1.4` deg | Sky field of view |
| `npix_x`, `npix_y` | `256` | Spatial pixels per axis |
| `map_output_file` | `'./Lco_cube_trial'` | Output path (no extension; `.npz` is appended) |
| `plot_cube` | `True` | Show map slice plot |
| `plot_pspec` | `True` | Show power spectrum plot |

The frequency bounds set the redshift range: `z_i = 115.27/34 − 1 ≈ 2.39`, `z_f = 115.27/26 − 1 ≈ 3.43` (COMAP low-z / "lowz" band).

The default catalogue path refers to a file not included in the repo — see `catalogues/` for available data.

---

## Subdirectory References

| Subdirectory | Reference doc | Contents |
|---|---|---|
| `limlam_mocker/` | [`limlam_mocker/CLAUDE_llm_llm.md`](limlam_mocker/CLAUDE_llm_llm.md) | Core Python package: halo loading, luminosity assignment, map binning, power spectrum |
| `catalogues/` | [`catalogues/CLAUDE_catalogues.md`](catalogues/CLAUDE_catalogues.md) | `.npz` catalogue format, `split_halo_catalogue.py`, COMAP mid-z metadata |
| `tables/` | [`tables/CLAUDE_tables.md`](tables/CLAUDE_tables.md) | `sfr_behroozi_release.dat` (SFR grid), `Tolgay2026_FIRE_CO.txt` (FIRE CO fits) |

---

## Authorship and Provenance

Original authors: George Stein (CITA) and Dongwoo Chung (Stanford). The inner package (`limlam_mocker/limlam_mocker/`) has been modified from the upstream version to:
- Accept `model_type` / `model_name` / `model_par` arguments routing to the `mass_luminosity.py` and `luminosity_functions.py` model libraries in the WebSkyLIM repo root (`Mhalo_to_Lline`).
- Support redshift limits up to z = 10 (extended from the original z = 4 cap) for the George Stein high-z COMAP data.
- Support `.h5` catalogues from Martine Lokken (PeakPatch centrals + satellites).

The upstream repository is `https://github.com/dongwooc/imapper2` (more feature-complete version by Dongwoo Chung and Tony Li).

---

## `.gitignore`

```
__pycache__
*/__pycache__
*.pyc
*.h5
*.npz
```

All compiled Python and catalogue data files are excluded. Catalogue `.npz` and `.h5` files must be obtained separately and placed in `catalogues/` or at whatever path is set in `params.py` / `params.catalogue_file`.
