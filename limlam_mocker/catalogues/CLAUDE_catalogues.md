# `limlam_mocker/catalogues/` — Reference

This directory holds halo-catalogue data files and the utility script for pre-processing them. The `.npz` and `.h5` catalogue files themselves are git-ignored (see `.gitignore`) and must be obtained separately; only metadata and tooling are committed here.

---

## `.npz` Catalogue Format

PeakPatch-output catalogues are stored as NumPy `.npz` archives. Each archive contains the following arrays:

| Key | Units | Description |
|---|---|---|
| `M` | M_sun | Halo mass |
| `x`, `y`, `z` | comoving Mpc | Cartesian position (lightcone axis along +z) |
| `vx`, `vy`, `vz` | km/s | Peculiar velocity |
| `zhalo` | — | Observed redshift (includes peculiar velocity contribution) |
| `zform` | — | Formation redshift of the halo |
| `cosmo_header` | — | Python dict (stored as a 0-d object array) with cosmological parameters |

**`cosmo_header` dict keys**: `h`, `Omega_M`, `Omega_B`, `Omega_L`, `ns`, `sigma8`. When a catalogue has been split into sub-fields by `split_halo_catalogue.py`, the header also contains `cen_x_fov` and `cen_y_fov` (degrees) — the angular offset of the sub-field centre from the catalogue origin, used by `load_halos.py` to re-centre RA/Dec at zero.

**Reading**: `np.load(filepath)` returns an `NpzFile` object. Access individual arrays with `['M']`, etc. Access the cosmo dict with `['cosmo_header'][()]` (the `[()]` unpacks the 0-d object array).

---

## Files

### `split_halo_catalogue.py`
Standalone script that takes a large single-pointing PeakPatch catalogue and tiles it into a grid of smaller sub-field catalogues. Designed to break a ~9.5 × 9.5 deg catalogue into a regular grid of 1 × 1 deg sub-fields.

**Hard-coded inputs** (edit before running):
- `filein = 'COMAP_z2.39-3.44_1140Mpc_seed_13579.npz'` — source catalogue
- `fov_x = fov_y = 9.52` deg — full catalogue angular size
- `fov_x_subfield = fov_y_subfield = 1.0` deg — output tile size

**What it does:**
1. Loads the full catalogue with a local `load_peakpatch_catalogue` function (mirrors `load_halos.py` but without a map instance dependency).
2. Computes `n_subfield = (fov_x // fov_x_subfield) × (fov_y // fov_y_subfield)` — number of output tiles (81 for the default 9×9 grid).
3. Loops over a 2D grid of tiles; for each tile calls `cull_peakpatch_catalogue` to select halos whose RA and Dec fall in `[fov_x_l, fov_x_r] × [fov_y_l, fov_y_r]`.
4. Saves each tile as `{original_name}_subfield_{i}.npz`, embedding the tile centre coordinates as `cen_x_fov` / `cen_y_fov` in `cosmo_header`.

**Output filenames**: `COMAP_z2.39-3.44_1140Mpc_seed_13579_subfield_0.npz` through `…_subfield_80.npz`.

This script executes on import (no `if __name__ == '__main__'` guard) — run it as `python split_halo_catalogue.py` from the `catalogues/` directory with the source `.npz` file present.

---

### `COMAP_z4.77-6.21_700Mpc/simdata_info.txt`
Metadata for the high-redshift COMAP simulation dataset (George Stein, CITA). Key facts:

- **Simulation name**: "midz"; frequency range 20–16 GHz (CO(1-0) from z ≈ 4.77–6.21)
- **Box size**: 700 Mpc comoving; halos span χ ≈ 7764–8464 cMpc
- **Geometry**: 18 independent realisations, each covering 4.72 × 4.72 deg; 72 non-overlapping ~4 sq deg patches total
- **Cosmology**: Li et al. parameters
- **Minimum halo mass**: 4×10⁹ M_sun (not M_sun/h — no h correction needed)
- **Format**: `.npz` compatible with `load_peakpatch_catalogue`
- **Original source**: `https://www.cita.utoronto.ca/~gstein/data/CO/COMAP_fullvolume/midz/z4.77-6.21_700Mpc/npz/`

The actual `.npz` files are not committed (git-ignored). To use them, download from the URL above and place in this subdirectory (or any path, then set `catalogue_file` in `params.py` accordingly).

---

## Notes for Navigation

- `load_halos.load_peakpatch_catalogue` (the package version) is what actually reads these files at runtime; `split_halo_catalogue.py` contains its own standalone copy of that loader.
- The `cen_x_fov` / `cen_y_fov` fields in sub-field headers are critical — `load_halos.py` subtracts them from the computed RA/Dec to re-centre the sub-field at (0, 0). Missing these causes the sky cut in `cull_peakpatch_catalogue` to silently discard all halos.
- Catalogue files are large (hundreds of MB per realisation for COMAP mid-z). The `.gitignore` excludes `*.npz` and `*.h5`.
