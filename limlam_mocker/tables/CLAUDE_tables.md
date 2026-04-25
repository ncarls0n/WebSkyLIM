# `limlam_mocker/tables/` — Reference

This directory holds tabulated data files used by the luminosity-assignment step of the simulation pipeline. Both files are read at runtime by functions in `halos_to_luminosity.py` (inner package) or `mass_luminosity.py` (repo root).

---

## `sfr_behroozi_release.dat`

**Source**: Behroozi, Wechsler & Conroy 2013a,b (arXiv:1207.6232, 1209.3013) — empirical star-formation rate as a function of halo mass and redshift, inferred by abundance matching.

**Grid dimensions**: 122 mass bins × 137 redshift bins = 16,714 rows (plus 2 header lines).

**Columns**:

| Col | Name | Units | Range |
|---|---|---|---|
| 1 | 1+z | — | 1.000 to 8.954 (z = 0 to 7.95) |
| 2 | log₁₀(M_halo) | M_sun | 9.0 to 15.05 |
| 3 | log₁₀(⟨SFR⟩) | M_sun/yr | −1000 (sentinel) to +2.19 |
| 4 | log₁₀(M_stellar) | M_sun | 5.84 to 10.75 |

**Sentinel value**: Cells where log₁₀(SFR) = −1000 represent physically inaccessible (mass, redshift) combinations — either below the resolution limit of the simulation or at very high redshift and very high mass where no halos exist in the Bolshoi/BolshoiP simulations used for the fitting. There are 3,710 such cells (~22% of the grid). `get_sfr_table(bad_extrapolation=True)` replaces these with a degree-4 spline extrapolation.

**How it is loaded** (`halos_to_luminosity.get_sfr_table`):
1. Read all four columns with `np.loadtxt`.
2. Convert: `dat_logzp1 = log₁₀(1+z)`, `dat_sfr = 10^log₁₀(SFR)` (i.e., linear SFR in M_sun/yr).
3. Reshape `dat_sfr` to `(122, 137)` — (mass axis, redshift axis).
4. Build a bilinear `scipy.interpolate.RectBivariateSpline` on `(log₁₀ M, log₁₀(1+z))`.
5. Query: `sfr_interp_tab.ev(log10(halos.M), log10(halos.redshift+1))` for vectorised evaluation.

**Used by**: `Mhalo_to_Lco_Li` (TonyLi 2016 model, legacy) — M → SFR → L_IR → L_CO. Not used by the modern `Mhalo_to_Lline` dispatcher, which routes to `mass_luminosity.py` models instead.

---

## `Tolgay2026_FIRE_CO.txt`

**Source**: Tolgay et al. 2026 (in prep.) — CO luminosity double-power-law fits calibrated against FIRE-2 zoom simulations. Produced by Tolgay Karademir. This is the data file backing the active-development `Tolgay2026CO` model in `mass_luminosity.py`.

**Format**: comma-delimited, 4 rows (one per redshift snapshot), no header.

**Columns**:

| Col | Symbol | Description |
|---|---|---|
| 1 | z | Redshift snapshot (0, 1, 2, 3) |
| 2 | A | Low-mass power-law slope (log-space) |
| 3 | B | High-mass power-law slope (log-space) |
| 4 | log C | log₁₀ of the normalisation at the characteristic mass M* |
| 5 | log M* | log₁₀(M*/M_sun) — the break/peak mass |
| 6 | σ | Log-normal scatter in dex |

**Full table**:

```
z,    A,      B,      log C,  log M*,  σ
0.0, -3.731, -1.249,  7.178,  11.037, 0.74
1.0, -3.357, -1.125,  7.957,  11.643, 0.80
2.0, -3.442, -2.119,  7.475,  11.558, 0.67
3.0, -3.175, -3.175,  7.205,  11.727, 0.54
```

**Physical interpretation of the double-power-law**:
```
log L_CO(M, z) = log C(z) + A(z) × log(M/M*(z))   for M < M*(z)
               = log C(z) + B(z) × log(M/M*(z))   for M > M*(z)
```
- A < 0: luminosity rises steeply with mass on the low-mass side (low-mass halos are inefficient CO emitters).
- B < 0, |B| < |A|: luminosity continues to rise but more shallowly at high mass (AGN feedback / CO excitation saturation).
- σ is the redshift-dependent intrinsic scatter applied as log-normal in `add_variable_log_normal_scatter`.

**How it is loaded** (`mass_luminosity.Tolgay2026CO`):
1. `np.genfromtxt(filepath, delimiter=',')` → shape (4, 6) array.
2. Parameters at each halo's redshift are obtained by linear interpolation along the z axis using `np.interp`.
3. Scatter is applied per-halo using `add_variable_log_normal_scatter` (which also reads this same file, identified by `interp_table[:,5]`).

**Note**: `add_variable_log_normal_scatter` in `halos_to_luminosity.py` contains a hardcoded absolute path `/home/njcarlson/Repositories/CITA_LIM/limlam_mocker/tables/FIRE-forged_CO.txt` as the default for this file. Always pass `interp_table_file` explicitly outside the original HPC environment, or use the `Tolgay2026CO` entry in `mass_luminosity.py` which resolves the path relative to the repo.

---

## Notes for Navigation

- Both files are read lazily at first use and cached in module-level globals (`sfr_interp_tab`) or local variables — they are not loaded at import time.
- The `sfr_behroozi_release.dat` path is resolved relative to the file location of `halos_to_luminosity.py` using `os.path.dirname(os.path.realpath(__file__))`, so it works regardless of the working directory.
- The `Tolgay2026_FIRE_CO.txt` path is similarly resolved in `mass_luminosity.py` (repo root), but the legacy `add_variable_log_normal_scatter` in `halos_to_luminosity.py` uses a hardcoded path — see warning above.
