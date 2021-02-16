# Spherical harmonic transforms in MATLAB

This repository provides a unified MATLAB interface to various spherical harmonic transform libraries.

`sht_plan(deg, n, backend)` creates a plan for a spherical harmonic transform using all spherical harmonics up to degree `deg` and an `n` x `n` grid in latitude and longitude. Transforms between spherical harmonic coefficients and values on the grid are done by the backend specified by `backend`. Currently, `backend` may be:

* `'fasttransforms'` - Mikael Slevinsky's FastTransforms library. (https://github.com/MikaelSlevinsky/FastTransforms)
* `'fmm3d'` - The Flatiron Institute's FMM3D library. (https://github.com/flatironinstitute/FMM3D)
* `'shtns'` - Nathanael Schaeffer's SHTns library. (https://bitbucket.org/nschaeff/shtns)

Each backend uses a different set of grid points on the sphere. The exact grid points used by a plan are stored in `plan.grid.lat` (for latitudinal grid points) and `plan.grid.lon` (for longitudinal grid points).

`sht_plan([lmax mmax], [nlat nlon], backend)` uses spherical harmonics with degree 0 ≤ l ≤ `lmax` and order |m| ≤ `mmax` and an `nlat` x `nlon` grid in latitude and longitude.

## Example
```
deg = 20;
n = 200;
c = randn((deg+1)^2, 1);
plan = sht_plan(deg, n, 'shtns');
V = plan.coeffs2vals(c);
d = plan.vals2coeffs(V);
```
