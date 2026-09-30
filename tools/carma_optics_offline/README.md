# Offline CARMA optics files for the alumina models

This tool writes CARMA's RRTMG optics files (`<model>_CRALUMnn_rrtmg.nc`) without running CAM. It uses CARMA's own Mie code, compiled offline.

## Why

CARMA normally writes these files at startup (`CARMA_CreateOpticsFile` in `cam/carma_intr.F90`), using its default Mie routine `miess` (Toon & Ackerman 1981). `miess` caps its series at |m|x terms, with a minimum of 135. Alumina has |m| < 1 in LW06 and LW07, next to its reststrahlen band, so `miess` stops with "The upper limit for acap is not enough" for radii above about 175 um.

The `alumina_d<D>um` particle size models therefore run with `carma_do_optics = .false.`. build-namelist points their `rad_climate` and `rad_diag` entries at `$DIN_LOC_ROOT/atm/cam/physprops/alumina/`.

## What it does

- `build.sh` compiles the pieces in this CAM tree into `carma_optics_driver`:
  - CARMA's `base/` code
  - CAM's `cam/carma_precision_mod.F90` and `cam/carma_constants_mod.F90`
  - CESM's `shr_kind_mod` and `shr_const_mod`
  - `cam_stubs.F90`, which stands in for `physconst`, `radconstants` and `cam_history_support`, using CESM's values
- `carma_optics_driver.F90` repeats what `carma_register` and `CARMA_CreateOpticsFile` do. It computes the band centres, creates the group, element and bins, then calls `getwetr` and `mie()` with the same unit conversions.
- `make_alumina_optics.py` reads `NBIN`, `rmin`, `vmrat`, the density, `mie_rh` and the `refidx` table from each model's `carma_model_mod.F90`. It then runs the driver and writes the files in `CARMA_CreateOpticsFile`'s layout, with provenance attributes added.

## Mie routine

Files are made with CARMA's `bhmie` (Bohren & Huffman), with one fix in `bhmie_patched.F90`. The downward recurrence for d(n) starts at `2*max(xstop,|mx|)+15` instead of `max(xstop,|mx|)+15`. The original start is not converged for large absorbing particles and gives up to 0.3% error in Qsca. The CAM source is not changed. Use `--mie toon` to run `miess` unchanged, for validation.

## Validation (2026-09-30)

- `--mie toon` on the base `alumina` model is bit-for-bit identical to the 36 files CARMA wrote at runtime.
- Patched `bhmie` against `miess`, wherever `miess` works:
  - base alumina, 36 bins: <= 1.2e-6
  - d2, d6 and d10um, all bins: <= 5e-5
  - d250um bins 1-2: <= 7e-6
- `bhmie` results do not change when the recurrence start is raised further (4x, or n = 10000).

## Usage

```
./build.sh
./make_alumina_optics.py --all-sizes            # the 10 alumina_d<D>um models
./make_alumina_optics.py alumina_d50um --outdir /path/to/inputdata/atm/cam/physprops/alumina
```

Rerun it whenever a model's bins, density or refractive indices change.
