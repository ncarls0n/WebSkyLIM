# `sfr_tables/` — Reference

This directory contains tabulated star-formation rate (SFR) data used by `mass_luminosity.py` models that require SFR(M, z) as an intermediate step. Both files are loaded from the repo root working directory — `mass_luminosity.py` references them by filename alone (not by path relative to `__file__`), so they must be present in whichever directory Python is invoked from.

---

## `sfr_release.dat` — Behroozi et al. 2013 SFR Table

**Source**: Behroozi, Wechsler & Conroy 2013a,b (arXiv:1207.6232, 1209.3013). Available from peterbehroozi.com/data.

**Content**: Identical to `limlam_mocker/tables/sfr_behroozi_release.dat` except this copy has no header lines — the data begins on line 1. The `limlam_mocker` copy has two comment lines at the top.

**Format**: 16,714 rows, 4 whitespace-delimited columns, no header.

| Col | Quantity | Units | Range |
|---|---|---|---|
| 1 | 1 + z | — | 1.000 to 8.954 (z = 0 to 7.95) |
| 2 | log₁₀(M_halo) | M_sun | 9.0 to 15.05 |
| 3 | log₁₀(⟨SFR⟩) | M_sun/yr | −1000 (sentinel) to +2.19 |
| 4 | log₁₀(M_stellar) | M_sun | 5.84 to 10.75 |

**Grid**: 122 unique mass bins × 137 unique redshift bins. The sentinel value −1000 marks physically inaccessible (mass, z) combinations (low-resolution or absent halos in the Bolshoi/BolshoiP simulations); 3,710 of 16,714 cells carry this value.

**How it is loaded** in `mass_luminosity.TonyLi` and `mass_luminosity.Behroozi_SFR`:
```python
x = np.loadtxt(BehrooziFile)          # BehrooziFile = 'sfr_release.dat'
zb     = np.unique(x[:,0]) - 1.       # redshift array (137 values)
logMb  = np.unique(x[:,1])            # log mass array (122 values)
logSFRb = x[:,2].reshape(137, 122, order='F')   # note: Fortran (column-major) reshape
logSFR_interp = interp2d(logMb, zb, logSFRb, bounds_error=False, fill_value=0.)
```
Note the `order='F'` reshape — the file is ordered with mass varying fastest (inner loop), so a Fortran-order reshape correctly produces the (z, M) grid. `bounds_error=False, fill_value=0.` means that out-of-range queries (very high z or very high M) return log(SFR) = 0, i.e., SFR = 1 M_sun/yr — a silent extrapolation that can cause errors for extreme halo masses.

**Used by**:
- `mass_luminosity.TonyLi` — `BehrooziFile` key in `MLpar` dict, default `'sfr_release.dat'`
- `mass_luminosity.Behroozi_SFR` — called directly by `SilvaCII` via `Behroozi_SFR('sfr_release.dat', Mvec, z)`
- `mass_luminosity.COMAP_Fid` (same TonyLi chain)

**Relationship to the `limlam_mocker` copy**: The two copies diverge only in the presence of a 2-line header. The `limlam_mocker` version (`sfr_behroozi_release.dat`) is loaded by `np.loadtxt(..., unpack=True)` which skips `#`-prefixed comment lines automatically. The `sfr_tables/sfr_release.dat` copy is loaded by `np.loadtxt` without comment handling — the lack of a header is therefore required for this loading code path. Do not add header lines to this file.

---

## `Silva15_SFR_params.dat` — Silva et al. 2015 SFR Parameters

**Source**: Silva, Vallini, Granato et al. 2015 (arXiv:1410.7315), Table 2. SFR(M, z) parameterised as a double power law in halo mass.

**Format**: 6 rows × 21 columns; whitespace-delimited; no header; scientific notation.

**Layout** — each row is one parameter, each column is one redshift snapshot:

| Row | Parameter | Units | Physical meaning |
|---|---|---|---|
| 1 | z | — | Redshift nodes: 0, 0.25, 0.5, 0.51, 1.07, 1.63, 2.19, 2.75, 3, 4, 5, 6, 7, 8.25, 10, 12, 13, 14.75, 16.5, 18.25, 20 |
| 2 | M₀ | M_sun/yr | SFR normalisation (values ~10⁻⁹ at z=0 to ~10⁻⁴ at z=20) |
| 3 | Mₐ | M_sun | Lower characteristic mass (all 10⁸ M_sun — constant in z) |
| 4 | M_b | M_sun | Upper characteristic mass (decreases from 8×10¹¹ at z=0 to ~1.6×10¹¹ at high z) |
| 5 | a | — | Low-mass slope (~2.7 at low z, ~2.1 at high z) |
| 6 | b | — | High-mass slope (negative; steepens from −4 at low z to ~−2.4 at high z, reflecting AGN/stellar feedback quenching) |

**Model equation** (Silva et al. 2015, Eq. 8):
```
SFR(M, z) = M₀(z) × (M/Mₐ(z))^a(z) × (1 + M/M_b(z))^b(z)
```
This is a double power law: rises as M^a for M ≪ Mₐ, transitions to M^(a+b) for Mₐ ≪ M ≪ M_b, then falls as M^(a+b) × (M/M_b)^b for M ≫ M_b. Since b < 0 and a > 0, the function has a peak and quenches at high mass.

Note: Mₐ = 10⁸ M_sun for all redshifts means the low-mass power law persists over many orders of magnitude in M. The high-mass break M_b drops from 8×10¹¹ to 1.6×10¹¹ M_sun between z=0 and z=6, reflecting earlier quenching in more massive halos at high redshift.

**How it is loaded** in `mass_luminosity.Silva_SFR`:
```python
x = np.loadtxt('Silva15_SFR_params.dat')
z0 = x[0,:]               # 21 redshift nodes
M0 = x[1,:]*u.Msun/u.yr   # normalisation
Ma = x[2,:]*u.Msun         # lower char. mass
Mb = x[3,:]*u.Msun         # upper char. mass
a  = x[4,:]               # low-mass slope
b  = x[5,:]               # high-mass slope
# Evaluate SFR at each z-node for all M, then interpolate in z:
for ii in range(z0.size):
    SFR0[ii,:] = M0[ii]*(M/Ma[ii])**a[ii]*(1+M/Mb[ii])**b[ii]
return interp1d(z0, SFR0.value, axis=0)(z)*u.Msun/u.yr
```
This evaluates the double power law at each of the 21 tabulated redshifts, stores the result as a (21, N_halo) grid, then uses `interp1d` to linearly interpolate to the actual halo redshifts.

**Used by**: `mass_luminosity.Silva_SFR` — which is called by `mass_luminosity.SilvaCII` in an older code path (currently commented out in favour of `Behroozi_SFR`; see the docstring of `SilvaCII`). The current active `SilvaCII` uses `Behroozi_SFR('sfr_release.dat', ...)` instead.

**Path dependency**: `Silva_SFR` loads with `np.loadtxt('Silva15_SFR_params.dat')` — a bare filename with no path resolution. The file must be present in the current working directory (i.e., the repo root) when `Silva_SFR` is called.

---

## Notes for Navigation

- `sfr_release.dat` is the live dependency for `TonyLi`, `COMAP_Fid`, and `SilvaCII` (via `Behroozi_SFR`). It must stay in the repo root working directory.
- `Silva15_SFR_params.dat` backs the `Silva_SFR` function which is currently **not called** by any active model — `SilvaCII` was updated to use `Behroozi_SFR` instead. The file is retained for the dormant code path.
- The `sfr_release.dat` here and `limlam_mocker/tables/sfr_behroozi_release.dat` carry the same data; only the header differs. If you update the Behroozi table (e.g., to Behroozi et al. 2019), update **both** copies.
- `sfr_reinterp.dat` (at the repo root, not in this directory) is a third copy re-interpolated onto a finer grid for the `LichenCII` model — it is referenced as `BehrooziFile` in the `Lichen`/`Lichen_v4` parameter presets in `params.py`.
