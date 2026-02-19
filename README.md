# WebSky Line Intensity Mapping

This repository contains the latest, stable version of the line intensity mapping (LIM) implementation of WebSky. These codes are based on the `limlam_mocker` package (originally written by George Stein, and later modified significantly by Dongwoo Chung, Patrick Breysse, Clara Chung and Patrick Horlaville). This implementation was forked from my branch of Patrick Horlaville's `CITA_LIM` repo and streamlined to include just critical components. Currently, this code is being developed by Nate Carlson and Doga Tolgay to produce FIRE-informed CO mocks.

## Summary of code components

The primary interface for interacting with the package is `lim.py` is the primary interface for interacting with the package. A lim model is loaded in python by importing `from lim import lim` and creating an object, e.g. `model = lim( model_params='COMAP_Fid', doSim=True )`

WebSkyLIM provides a lightweight LIM toolkit with three tightly-coupled layers:

1. **Astrophysical modeling (`LineModel`)**  
   Compute luminosity functions or mass–luminosity models, mean intensity, bias, and analytic LIM power spectra.

2. **Instrument + survey forecasting (`LineObs`)**  
   Add survey geometry and instrument noise to predict quantities like voxel noise, noise power, k-space errors, and SNR.

3. **Map-based simulations (`LimLam`)**  
   Wrap the `limlam_mocker` pipeline to generate 3D LIM cubes from halo catalogues and measure map-based statistics (e.g., simulated power spectra, VID/histograms). These are built from halo catalogues using response functions in `mass_luminosity.py`.

Everything is built around `astropy.units`, so **inputs should carry units** and outputs come back with units whenever possible.

---

## Repository structure

### Core interface
- `lim.py`  
  Primary user-facing interface. Provides the factory function:
  ```python
  from lim import lim
  m = lim(model_params='COMAP_Fid', doSim=True)
  ```

### Core model/forecast classes
- `line_model.py`  
  Defines `LineModel`: cosmology + line model (LF or ML) + analytic statistics.
- `line_obs.py`  
  Defines `LineObs(LineModel)`: adds instrument/survey parameters and forecasting utilities.
- `limlam.py`  
  Defines `LimLam(LineObs)`: runs/reads `limlam_mocker` simulations and exposes map-based stats.

### Model libraries
- `mass_luminosity.py`  
  Mass→luminosity response functions for “ML” models (e.g., CO/CII models, FIRE-informed options, etc.).
- `luminosity_functions.py`  
  LF parameterizations for “LF” models (e.g., Schechter-like forms, file-based LF interpolation).

### Parameter presets / configs
- `params.py`  
  A collection of ready-to-use parameter dictionaries (instrument + cosmology + model choices).  
  **Important:** some presets include machine-specific absolute paths for halo catalogues — you will need to update `catalogue_file` for your environment.

### Data tables / helpers
- `sfr_tables/`, `sfr_release.dat`, `sfr_reinterp.dat`  
  Star-formation history tables used by some models.
- `fire_galaxies/`  
  FIRE-related files used in FIRE-informed modeling workflows.
- `_utils.py`, `_vid_tools.py`  
  Shared utilities + VID/CVID helper functions.
- `stacking*.py`, `stacking_utils.py`, `stacking_params.py`  
  Stacking pipelines/utilities (work-in-progress/analysis scripts).
- `limlam_mocker/`  
  Included `limlam_mocker` codebase used for map simulations.

---

## Installation

There is currently no pinned `requirements.txt` in this repo, so install dependencies manually.

### Option A: pip + venv (recommended)
```bash
git clone https://github.com/ncarls0n/WebSkyLIM.git
cd WebSkyLIM

python -m venv .venv
source .venv/bin/activate

pip install --upgrade pip
pip install numpy scipy astropy hmf h5py multiprocess
```

### Make the repo importable
From the repo root, either:
- run scripts from this directory, **or**
- add the repo to `PYTHONPATH`, e.g.
  ```bash
  export PYTHONPATH="$(pwd):$PYTHONPATH"
  ```

> Notes:
> - `LimLam` can also read PeakPatch-formatted outputs. If you use that path, you may need `peakpatchtools` available in your environment.
> - Some workflows may assume an HPC environment and/or large halo catalogues.

---

## Quickstart

### 1) Analytic model only (no instrument)
```python
import astropy.units as u
from lim import lim

m = lim(model_params='default_par', doObs=False)

print("z =", m.z)
print("Tmean =", m.Tmean)
print("Pshot =", m.Pshot)
print("Pk[0:5] =", m.Pk[:5])
```

### 2) Add instrument + survey forecasting (recommended default)
```python
from lim import lim

m = lim(model_params='default_par', doObs=True, doSim=False)

print("Nvox =", m.Nvox)
print("sigma_N =", m.sigma_N)
print("Pnoise =", m.Pnoise)
print("SNR =", m.SNR)
```

### 3) Run a simulation-backed cube (requires a halo catalogue)
```python
import astropy.units as u
from lim import lim

# Start from a preset, but override the catalogue path
m = lim(
    model_params='COMAP_Fid',
    doObs=True,
    doSim=True,
)

# If needed, override catalogue location after construction
m.update(catalogue_file="/path/to/your/catalogue.npz")

cube = m.maps  # (Nx, Ny, Nch) with units
print("Cube shape:", cube.shape)
print("Simulated mean:", m.Tmean_sim)

# Simulated power spectrum
k = m.k_sim
Pk = m.Pk_sim
```

---

## The `lim()` factory

`lim.py` exposes a single function that returns one of the following objects depending on flags:

- `LineModel` if `doObs=False`
- `LineObs` if `doObs=True` and `doSim=False`
- `LimLam` if `doObs=True` and `doSim=True`

```python
from lim import lim
m = lim(model_params='TonyLi_PhI', doSim=True)
```

### Inputs
- `model_params`  
  Either:
  - a **string** naming a dictionary in `params.py` (e.g. `"COMAP_Fid"`), or
  - a **dict** containing parameter names/values directly.
- `doObs`  
  Include instrument/survey modeling.
- `doSim`  
  Enable `LimLam` simulation interface.
- `match_sim_cosmo`  
  If simulations and analytic cosmology differ, optionally update analytic cosmology to match the simulations.

---

## Choosing a line model

Two model types are supported:

### A) LF models (`model_type='LF'`)
Define an **emission luminosity function** `dn/dL`, and compute LIM quantities by integrating over luminosity.

- Add new LF models in `luminosity_functions.py`
- Select them via:
  ```python
  model_type = 'LF'
  model_name = 'SchCut'   # e.g. Schechter with low-L cutoff
  model_par  = {...}
  ```

### B) ML models (`model_type='ML'`)
Define **luminosity as a function of halo mass** `L(M)` (optionally with scatter/duty cycle), and compute LIM quantities by integrating over the halo mass function.

- Add new ML models in `mass_luminosity.py`
- Select them via:
  ```python
  model_type = 'ML'
  model_name = 'TonyLi'   # example
  model_par  = {...}
  ```

### Scatter and duty cycle
Many workflows use:
- `sigma_scatter` : log10 scatter in luminosity at fixed mass (mean-preserving)
- `fduty` : duty cycle fraction

These are available at the `LineModel`/`LineObs` level and also used in `LimLam` map-making.

---

## Working with presets (`params.py`)

`params.py` contains a mix of “generic defaults” and instrument/model presets (COMAP-like, CCATp-like, CII/CO variants, etc.). Examples include:
- `default_par` (baseline defaults)
- `TonyLi_PhI`, `TonyLi_PhII` (Tony Li style models)
- `COMAP_Fid` (COMAP fiducial-style preset)
- multiple `Silva_*` and `Lichen_*` presets (CII/CO variations)

**Important:** Some entries set `catalogue_file` to an absolute path used during development. You should override:
```python
m.update(catalogue_file="/your/path/catalogue.npz")
```

---

## Outputs you’ll commonly use

### From `LineModel` (analytic)
- `z` : target redshift from `nu / nuObs - 1`
- `Tmean` : mean brightness temperature (or Jy/sr if `do_Jysr=True`)
- `Pshot` : shot-noise term
- `Pk` : analytic 3D LIM power spectrum
- `bavg`, `nbar`, `dndM`, `dndL`, `LofM` (model-dependent)

### From `LineObs` (instrument + survey)
- Geometry:
  - `Nch`, `Nside`, `Npix`, `Nvox`
  - `Vfield`, `Vvox`
- Noise:
  - `sigma_N` (per voxel)
  - `Pnoise`
- Sensitivity:
  - `sk` (total uncertainty per k-bin)
  - `SNR`

Also includes support for:
- 1D power spectrum style quantities (`P1D`, `SNR_1D`, etc.)
- angular power spectrum predictions (`Cl`) in exact or Limber mode

### From `LimLam` (simulation/maps)
- `maps` : simulated cube (with units)
- `k_sim`, `Pk_sim`, `cv_sim` : simulated power spectrum + errors
- `Tmean_sim`, `Tmean_halos` : map mean and halo-sum mean
- `Pshot_halos` : halo-estimated shot noise
- `Bi_sim`, `PofN_sim` : map-based VID / number-count PDFs
- `noise_added_map` : example helper to add noise to a beam-smoothed cube

---

## Updating parameters safely

These classes use cached properties for speed. **Do not mutate attributes directly.**  
Always call `update()` so cached quantities are cleared and recomputed correctly.

```python
import astropy.units as u
m.update(
    nuObs = 30*u.GHz,
    beam_FWHM = 4.1*u.arcmin,
    sigma_scatter = 0.3,
)
```

---

## Halo catalogues and formats (for simulations)

`LimLam` supports:
- `.npz` catalogues (common for PeakPatch-derived products)
- `.h5` catalogues
- PeakPatch output formats (requires PeakPatch tools; see `limlam.py`)

At minimum you will need to provide:
- `catalogue_file` : location of the halo catalogue
- optionally `map_output_file` : where `limlam_mocker` writes outputs

---

## Adding new models

### Add a new LF model
1. Implement a function in `luminosity_functions.py` with signature:
   ```python
   def MyLF(Lvec, LFparams):
       ...
       return dndL
   ```
2. Use it by setting:
   ```python
   model_type = 'LF'
   model_name = 'MyLF'
   model_par  = {...}
   ```

### Add a new ML model
1. Implement a function in `mass_luminosity.py` with signature:
   ```python
   def MyML(M, model_par, z):
       ...
       return L
   ```
2. Use it by setting:
   ```python
   model_type = 'ML'
   model_name = 'MyML'
   model_par  = {...}
   ```

Tip: Add a new preset dict in `params.py` once your model works, so it’s easy to reproduce.

---

## Known gotchas / practical notes

- **Absolute paths in presets:** Several `params.py` presets were created on specific machines and include absolute paths (especially `catalogue_file`). Override them.
- **Units everywhere:** Inputs should carry `astropy.units`. This is intentional.
- **Cached properties:** Use `update()` instead of attribute mutation.
- **Simulation cosmology mismatch:** If your analytic cosmology doesn’t match a catalogue’s cosmology, the code can optionally set analytic cosmology to match simulations via `match_sim_cosmo=True`.

---

## Contributing

Issues and PRs are welcome. If you add:
- a new line model,
- a new instrument preset,
- new catalogue readers,

please include:
- a short description,
- a minimal example snippet,
- and (ideally) a `params.py` preset that reproduces your test.

---

## Acknowledgements

This repository builds on the LIM ecosystem developed across several projects and groups, particularly `limlam_mocker` and WebSky-related workflows, with major contributions and modifications by George Stein, Dongwoo Chung, Patrick Breysse, Clara Chung, and Patrick Horlaville.

---

## License

A license file is not currently included in the repository. If you intend this to be used broadly, consider adding a `LICENSE` file (e.g., MIT/BSD/GPL) and a citation guideline.
