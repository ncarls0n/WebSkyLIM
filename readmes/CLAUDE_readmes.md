# `readmes/` — Reference

This directory contains two Jupyter notebooks written by Patrick Horlaville (phorlaville24@ubishops.ca) that document the usage of the WebSkyLIM codebase. They serve as tutorials and worked examples. Both notebooks have hardcoded HPC paths and cannot be run locally without modification.

The companion paper for both notebooks: arXiv:2309.15733 (Horlaville et al. 2023).

---

## `llm_readme.ipynb` — Getting Started with `limlam_mocker`

**Purpose**: Step-by-step tutorial for generating line intensity maps using the `lim()` class. Covers the most important `lim` object properties and shows how to visualise and interpret the outputs.

**Audience**: New users; read before `dsrel_readme.ipynb`.

### What it demonstrates

**Setup**: Instantiates `lim('Lichen_v4', doSim=True)` then calls `m_cii.update(...)` to set all parameters explicitly. The `update()` call is the recommended pattern — it propagates parameter changes through the `cached_property` invalidation system.

**Parameters shown** (with example values):
| Parameter | Example | Description |
|---|---|---|
| `model_par` | `{'zdex':0.4, 'M0':1.9e9, 'Mmin':2e10, 'alpha_MH1':0.74, 'alpha_LCII':0.024, ...}` | `LichenCII` model parameters |
| `dnu` | `2.8 GHz` | Channel bandwidth |
| `nuObs` | `270 GHz` | Central observed frequency (CII at z ≈ 5.7) |
| `Delta_nu` | `40 GHz` | Total survey bandwidth |
| `tobs` | `40000 h` | Total observation time |
| `Omega_field` | `4 deg²` | Survey sky area |
| `catalogue_file` | HPC path to `.npz` | PeakPatch halo lightcone |

**Key outputs demonstrated**:

- `m_cii.maps` — The 3D intensity cube (shape: `(Nside, Nside, Nch)`), in Jy/sr. Loading takes ~2 min for a full lightcone.
- `m_cii.mapinst.nu_bincents` — Frequency channel centres; used to recover per-slice redshift via `z = nu_rest/nu - 1`.
- `m_cii.noise_added_map` — Signal + injected Gaussian noise cube (σ = `m_cii.sigma_N` per voxel).
- `m_cii.beam_width`, `m_cii.Nside`, `m_cii.Omega_field` — Beam and survey geometry for Gaussian filtering.

**Beam convolution pattern**:
```python
from scipy.ndimage import gaussian_filter
# 1-pixel beam (default resolution):
gaussian_filter(cii_cube[:,:,slice_idx], 1)
# Custom beam size: use the pix_res() helper defined in the notebook
```

**`pix_res()` helper** (defined in notebook, not in any module): converts a physical beam FWHM to pixel width for use with `gaussian_filter`. Formula: `beam_width × map_dim_pix / map_dim_deg`.

**HPC path dependency**: `catalogue_file` is set to `/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz` — must be overridden for local use. The `%cd` cell also hardcodes a CITA Niagara path.

---

## `dsrel_readme.ipynb` — Relative Entropy Analysis

**Purpose**: Documents the point-wise differential relative entropy (PRE) analysis applied to intensity map voxel distributions — a technique for quantifying how sensitively the voxel intensity PDF changes with model parameter perturbations. Read `llm_readme.ipynb` first.

**Reference**: arXiv:2309.15733, §2.2 and §4.3.

### Physics: Differential Relative Entropy

The **relative entropy** (Kullback–Leibler divergence) between two voxel distributions P and Q is:
```
D_KL(P || Q) = Σ_i P_i × log(P_i / Q_i)
```
The **point-wise relative entropy (PRE)** per intensity bin is the integrand `P_i × log(P_i / Q_i)`, summed over bins to give D_KL. To measure how D_KL changes with a model parameter λ, the notebook computes the finite-difference derivative:
```
ΔPRE/Δλ ≈ P_baseline × [log(P_λ+Δ / P_baseline) − log(P_λ−Δ / P_baseline)] / (2Δλ)
```
This tells you which intensity bins (i.e., which voxel brightness range) are most discriminating for constraining each parameter.

### Parameters varied

The notebook sweeps these `LichenCII` model parameters around their fiducial values:

| Parameter | Fiducial | Range swept | Physical meaning |
|---|---|---|---|
| `zdex` | 0.4 | 0.3, 0.4, 0.5 | Log-normal scatter in z |
| `alpha_LCII` | 0.024 | 0.020, 0.024, 0.028 | CII luminosity per unit HI mass |
| `alpha_MH1` | 0.74 | 0.68, 0.74, 0.80 | HI mass – halo mass power-law slope |
| `log(Mmin)` | 10.0 | 9.9, 10.0, 10.1 | Log of minimum halo mass |
| `alpha0` | −1.412 | −1.312, −1.412, −1.512 | Halo bias parameter |
| `gamma0` | 0.31 | 0.27, 0.31, 0.35 | Halo bias redshift evolution |

### Two noise regimes

The notebook loads pre-computed voxel histograms for two cases:
- **Signal only** (`sig/`): pure CII intensity cube with no noise
- **Stage 2** (`for/for_tobs40kh/`): signal + instrument noise (t_obs = 40,000 h CCAT-prime Stage 2)

The `.npy` files (e.g., `n_avg_zdex0-3.npy`) are ensemble-averaged histograms over multiple lightcone realisations, stored on HPC scratch at `/mnt/scratch-lustre/horlaville/nuObs270/`. These are not committed to the repo.

### What the notebook cannot run locally

- The `%cd` magic points to a CITA HPC path.
- All `np.load(...)` calls load pre-computed `.npy` histograms from HPC scratch — none of these files exist locally.
- The single-lightcone PRE computed inline (Cells 15–25) does work if a catalogue file and `m_cii` instance are available, but the ensemble results (Cells 29–57) require the HPC data.

---

## Notes for Navigation

- Neither notebook is imported or called by any pipeline code — they are documentation only.
- The `pix_res()` helper in `llm_readme.ipynb` would be a useful addition to `line_obs.py` or `_utils.py` if beam-convolved maps become a routine output.
- The `sfr_reinterp.dat` file referenced in `m_cii.update(model_par={'BehrooziFile':'sfr_reinterp.dat'})` is at the repo root (not in `sfr_tables/`); it is a re-interpolated version of the Behroozi table on a finer grid for the `LichenCII` model.
