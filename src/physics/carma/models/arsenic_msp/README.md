# CARMA model `arsenic_msp`: reentry arsenic on meteoric smoke

Meteoric smoke (as in `meteor_smoke`: group `DUST`, 28 bins) carries an arsenic
core-mass element `DUSTAS`. A lumped arsenic gas `ASGAS` (mass of As atoms) is
emitted at 85 km and taken up by the smoke at the gas-kinetic collision rate with
a Fuchs-Sutugin correction and sticking coefficient `carma_as_gamma`. Design:
`~/MIT/cesm2_port/ARSENIC_MSP_CARMA_PLAN.md`.

## Files

| File | Kind | Purpose |
|---|---|---|
| `carma_model_mod.F90` | model | groups, elements, gas, emissions, `CARMA_EmitGas` |
| `carma_model_flags_mod.F90` | model | `carma_model_nl` namelist |
| `arsenic_uptake_mod.F90` | new, CARMA side | uptake of `ASGAS` onto smoke (`arsenic_uptake`) |
| `newstate_calc.F90` | override of CARMA `base/` | calls `arsenic_uptake` each CARMA step |
| `carma_intr.F90` | override of CAM `carma/cam/` | registers, initializes, emits and outputs the gas |
| `*.patch` | | `diff -u <original> <override>` for each override |

## How the overrides take effect

`bld/configure` writes the CARMA directories into `Filepath` in this order:
`carma/models/<model>`, `carma/cam`, `carma/base`. A file in the model directory
therefore replaces the file of the same name in `cam/` or `base/`. The CIME
Makefile adds `-I` for every `Filepath` directory, so `#include "carma_globaer.h"`
(in `base/`) resolves from here. The precedent is `models/alumina/vertical.F90`
on the alumina branch.

These overrides are used only when CAM is configured with `-carma arsenic_msp`.

## Overrides

Every change in an override is between `!++AS` and `!--AS`. Each file starts with
a header naming the original.

### `newstate_calc.F90`

- **Original:** `src/physics/carma/base/newstate_calc.F90`, CARMA SVN tag
  `carma3_49_rel` (r91872). `base/` is not in git; it is an SVN checkout.
- **Changes:**
  1. `use arsenic_uptake_mod, only: as_do_uptake, arsenic_uptake`.
  2. Immediately after the `if (do_coag) call microslow(...)` block, and before
     `if (do_grow)`: if `as_do_uptake`, call `arsenic_uptake` for `iz = 1, NZ`.
- **Why:** the natural place for gas-to-particle transfer is `microfast`, but
  `microfast` runs only `if (do_grow)`, and `arsenic_msp` keeps growth off. This
  spot runs once per call of `newstate_calc`, after coagulation and before the
  (unused) growth loop, whatever `do_grow` is.
- `newstate_calc` is called twice per step (in-cloud and clear sky) when
  `do_incloud` is on. `carma_model_mod` refuses that setting at initialization and
  `arsenic_uptake` returns an error if it sees it.

### `carma_intr.F90`

- **Original:** `src/physics/carma/cam/carma_intr.F90` at CAM tag
  `cam_cesm2_1_rel_60` (a03b84b7), 2723 lines.
- **Changes:**
  1. `carma_register`: the original calls `cnst_get_ind(shortname, icnst4gas(igas))`,
     which ends the run if the gas is not already a CAM constituent. It now calls it
     with `abort=.false.` (returns -1 when not found) and, if not found,
     `cnst_add(shortname, wtmol, cpair, 0._r8, icnst4gas(igas), longname=name,
     mixtype=carma_mixtype)`. A module array `carma_added_gas(NGAS)` records which
     gases were added here.
  2. `carma_init`, gas loop: `addfld` for `<gas>EM` (kg/kg/s), and, for gases added
     in change 1 only, `addfld` for the constituent itself (`ASGAS`, kg/kg; also an
     `add_default` when `history_carma`). In CAM each package that adds a
     constituent must `addfld` it; otherwise `ASGAS` cannot be written to history.
  3. `carma_emission_tend`, after the element loop: for each gas, call
     `CARMA_EmitGas(carma, igas, icnst, dt, state, cam_in, tendency, rc)` from
     `carma_model_mod`, put `tendency` into `ptend%q(:ncol,:pver,icnst)` and write
     `<gas>EM`.
  4. `carma_init_cnst`: if `name` is a gas constituent, set `q = 0` where `mask` is
     true. The original only handled elements, so `q` (`intent(out)`) was left
     undefined for gases.
  5. `carma_init` and `carma_timestep_tend`: history field `ASGASK` (1/s, `lev`,
     averaging flag 'A'), the arsenic gas loss rate `k` that `arsenic_uptake` used at
     each level (added 2026-10-21, user-approved). `carma_init` registers it after
     the gas loop. In `carma_timestep_tend`, `arsenic_uptake_reset(pver)` clears the
     module array before each column's `CARMASTATE_Step`, `arsenic_uptake_getk`
     copies it into a local `asgask(pcols,pver)` right after the step, and
     `outfld('ASGASK', ...)` writes the chunk after the column loop. Level mapping:
     the column is passed to `CARMASTATE_Create` as `pmid(icol,1:pver)` (CAM order,
     k = 1 at the top) and CARMA keeps that order for `I_HYBRID`; the interface maps
     results back with `ptend%q(icol,1:pver,:) = newstate(1:pver)`. So CARMA `iz`
     is CAM `k` and no flip is needed. Every column `1..ncol` goes through
     `CARMASTATE_Step` (none is skipped), so every written value is the column's own
     `k`. `ASGASK` is always registered and written; it is 0 where no rate was
     computed (`as_gamma <= 0` or `gc <= 0` at that level), and everywhere when
     `carma_as_do_uptake = .false.`. Where `k*dtime < 1e-30` the computed `k` is
     written although nothing is moved. No fill value is used.
  6. `carma_register`, right after `CARMA_Initialize`: if `carma_as_top_noflux`
     (default `.true.`, added 2026-10-24), set `carma%f_ibbnd_pc = I_FLUX_SPEC`.
     See "Model-top boundary" below.
  7. Bin-summed surface deposition (added 2026-10-24): for each element of a
     prognostic group, `carma_init` registers `<element>SW` and, if the group does dry
     deposition, `<element>DD` (kg/m2/s, horizontal, averaging flag 'A'). Here that is
     `DUSTASSW`, `DUSTASDD` (arsenic) and `DUSTSW`, `DUSTDD` (smoke plus arsenic, as
     `DUSTnn`). `carma_output_diagnostics` writes `<element>DD` as the sum of
     `sbdiags(:, 1:maxbin, ielem, SBDIAGS_DD)`; `carma_wetdep_tend` writes
     `<element>SW` as the sum of the per-bin `sflx`. Both are the sums of the per-bin
     fields over bins 1..maxbin, so production runs need 2 fields per element
     instead of 2 x 28. Neither is on a tape by default; add them to `fincl`.
- **Why:** the stock interface assumes CARMA gases (H2O, H2SO4) already exist as CAM
  constituents and have no CARMA emission. `ASGAS` is new and emitted by CARMA.
  `ASGASK` gives the exact `k` for WP7 test 2; an offline recomputation from the
  `:I` snapshot cannot see the state the uptake sees (see
  `Postprocessing/arsenic/README.md`).

## Uptake physics (`arsenic_uptake_mod.F90`)

For each level, in cgs:

- `n_air = rhoa/(xmet*ymet*zmet) * AVG / WTMOL_AIR` (as in `base/setupgrow.F90`)
- `D = 1.76575e17 * sqrt(T) / n_air`, `vbar = sqrt(8 RGAS T / (PI wtmol))`,
  `lambda = 3 D / vbar`
- per bin: `N = pc(iz,ibin,ienconc)/(xmet*ymet*zmet)`, `r = r_wet`, `Kn = lambda/r`,
  `f = (1+Kn)/(1 + (4/(3 gamma)+0.377) Kn + 4/(3 gamma) Kn^2)`,
  `k_i = 4 PI r D N f` (or `gamma PI r^2 vbar N` with `carma_as_do_fuchs=.false.`)
- `dgc = gc (1 - exp(-k dtime))` is moved from `gc` to `pc(:,:,iecore)` in
  proportion to `k_i / k`. A series is used for `1 - exp(-x)` when `x < 1e-4`.
- Returns early if `gamma <= 0`, `gc <= 0` or `k dtime < 1e-30`.
- Diagnostic: the `k` used at each level is stored in the module array `as_kused`
  (`!$OMP THREADPRIVATE`, since CAM may run chunks on several threads). Each call
  first sets `as_kused(iz) = 0` and then, once the sum over bins is done, `ksum`
  (also when `k dtime < 1e-30`). `arsenic_uptake_reset(nz)` (zero, allocate if
  needed) and `arsenic_uptake_getk(nz, k)` are called from the `carma_intr.F90`
  override; it becomes history field `ASGASK`. The diagnostic does not change
  any prognostic value.

`gc` and a core-mass `pc` are both mass mixing ratio times `rhoa_wet`
(`CARMASTATE_SetGas`/`GetGas`, `CARMASTATE_SetBin`/`GetBin` in
`base/carmastate_mod.F90`), so the transfer needs no unit conversion and conserves
mass to round-off.

### Bin shift (`carma_as_do_binshift`, default `.true.`, added 2026-10-22)

The arsenic taken up adds to the particle mass (risk R2 in the plan; the user chose
this over a passive label on 2026-10-22). CARMA gives every particle in bin `i` the
mass `rmass(i)` (smoke bins double in volume, `vmrat = 2`). After the uptake, in the
same call and level:

- `m = rmass(i) + dA_i / N_i`, where `dA_i` is the arsenic just added to bin `i`
  and `N_i` is the number element (the metric scaling cancels).
- `j` = the largest bin with `rmass(j) <= m`. A fraction
  `fmove = (m - rmass(j)) / (rmass(j+1) - rmass(j))` of the particles goes to bin
  `j+1` and the rest to bin `j`. When `j = i` (the usual case) `fmove` is computed
  from `dA_i / N_i` directly, to avoid cancellation.
- The arsenic core of bin `i` (old core plus `dA_i`) is split with the same
  fractions. All bins are moved from the state before the shift.
- Number, total particle mass (`N rmass`) and arsenic mass are conserved to
  round-off. Particles that would grow past bin `NBIN` (radius about 0.1 um) stay
  there, losing the mass above `rmass(NBIN)`; that bin holds almost no smoke.
- This is first-order upwind in mass space, so it adds some numerical spreading
  across bins. Each step moves only a small fraction of a bin.

Consequences:

- `DUSTnn` is smoke **plus** arsenic. Smoke-only mass is `DUSTnn - DUSTASnn`.
- Core mass can no longer exceed the particle mass because of the uptake, so the base
  `rhopart.F90` repair branch acts only on advection noise and on the
  `smallconc` core floor.
- With `carma_do_fixedinit = .true.` the fall speeds and coagulation kernels stay
  those of bare smoke for each bin. The shift captures the larger size (mass), not
  the higher density of arsenic (5.73 vs 2.0 g/cm3). `carma_do_partialinit` would
  add the density effect on fall speeds (not tested).
- With `carma_as_do_binshift = .false.` the model is the v1 passive label: `DUSTnn`
  is smoke only and particles never change bin because of arsenic.

Science notes: `D` uses the kinetic-theory formula CARMA uses for H2SO4; the
uncertainty this adds is probably small next to the `gamma = 1` assumption. The
module variables (`as_do_uptake`, `as_do_fuchs`, `as_do_binshift`, `as_gamma`, `as_igas`,
`as_igroup`, `as_iecore`) are set by `CARMA_InitializeModel`.

Known limitation: CARMA treats a core as part of the bin mass `rmass`, not as
extra mass. `base/rhopart.F90` computes the particle density from
`rmass*N - mcore` (smoke) plus `mcore` (arsenic, 5.73 g/cm3), so arsenic raises
`rhop` and the fall speed of the smoke. Where core mass exceeds `rmass*N` it caps
the density and resets the number, `pc(iepart) = mcore/rmass`. The uptake therefore
can change `DUSTnn`; see risk R2 in the plan.

## Model-top boundary (`carma_as_top_noflux`, default `.true.`, added 2026-10-24)

`CARMA_InitializeVertical` (`base/carma_mod.F90`) hardcodes both particle-transport
boundaries to `I_FIXED_CONC` with a concentration of 0. For `I_HYBRID`, level 1 is the
model top, and `vertadv`/`vertdif` apply `ibbnd_pc` at the level-1 interface (and
`itbnd_pc` at `NZP1`, the surface). So molecular and eddy diffusion carried smoke and
arsenic out through the model top into a zero-concentration boundary.

`base/vertical.F90` sets `pc_surf` (the `DD` history fields) to the column loss in each
transport step, which is only the surface flux if nothing leaves through the top. With
the default boundary, `DD` was dominated by this top loss. In the WP8 2-day run,
`DUSTASnnDD` was about 1e-23 kg/m2/s, almost all in bin 1, while `DUSTAS` in the lowest
level was about 1e-48 kg/kg.

With the flag on, the level-1 boundary is `I_FLUX_SPEC` with `fbotpart = 0`: `vertadv`
sets `vtrans(1) = 0` and `vertdif` sets `vertdifu(1) = vertdifd(1) = 0`, so there is no
flux through the top and `DD` is the surface flux only. The surface boundary is
unchanged. This also affects the smoke, which is no longer bit-for-bit with the CARMA
default. `.false.` restores the CARMA default.

The change is in this model's `carma_intr.F90` override, not `CARMA_InitializeModel`,
because `CARMA_InitializeModel` gets `carma` as `intent(in)`, and `CARMA_Initialize`
(which calls `CARMA_InitializeVertical`) runs after `CARMA_DefineModel`.

## Regenerating the patches

From `src/physics/carma` in the CAM checkout (`base/` must be the
`carma3_49_rel` SVN checkout):

```sh
M=models/arsenic_msp
diff -u --label base/newstate_calc.F90 --label $M/newstate_calc.F90 \
    base/newstate_calc.F90 $M/newstate_calc.F90 > $M/newstate_calc.F90.patch
diff -u --label cam/carma_intr.F90 --label $M/carma_intr.F90 \
    cam/carma_intr.F90 $M/carma_intr.F90 > $M/carma_intr.F90.patch
```

`diff` exits with status 1 when the files differ; that is expected.

## Rebasing onto a new CAM or CARMA version

1. Copy the new original over the override, for example
   `cp cam/carma_intr.F90 $M/carma_intr.F90`.
2. Apply the old patch: `patch $M/carma_intr.F90 < $M/carma_intr.F90.patch`.
   Resolve any rejected hunks by hand, keeping every `!++AS`/`!--AS` block, and
   update the version named in the header comment.
3. Check that the upstream changes do not make an override unnecessary (for
   example, if upstream `carma_intr` starts adding gases with `cnst_add`).
4. Regenerate the patch as above, rebuild, and rerun the WP7 tests (at least mass
   closure and the smoke-unchanged test).
