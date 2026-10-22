# The CARMA alumina model in CESM2.1.5 / WACCM: overview and how-to

This guide is for someone new to CARMA and CESM. It covers:

1. [Background: CESM, CAM/WACCM and CARMA](#1-background-cesm-camwaccm-and-carma)
2. [How the alumina model works](#2-how-the-alumina-model-works)
3. [What was changed in CAM to add it](#3-what-was-changed-in-cam-to-add-it)
4. [Installation: getting the code](#4-installation-getting-the-code)
5. [Input files](#5-input-files)
6. [Running a case, step by step](#6-running-a-case-step-by-step)
7. [`user_nl_cam` reference](#7-user_nl_cam-reference)
8. [Output fields and analysis tips](#8-output-fields-and-analysis-tips)
9. [Adding a new particle-size model](#9-adding-a-new-particle-size-model)
10. [Troubleshooting](#10-troubleshooting)

The code is on the branch **`ashajain/carma_aluminum_model`** of **https://github.com/ashajain652/CAM_carma19** (a fork of NCAR's CAM, based on `cam_cesm2_1_rel`). The scripts that make the emission and size-distribution input files are in **https://github.com/ashajain652/alumina_emissions**.

---

## 1. Background: CESM, CAM/WACCM and CARMA

- **CESM** is a coupled Earth system model. Each part (atmosphere, land, ocean, sea ice and so on) is a "component" with its own git repository. The top-level CESM checkout lists the components and their versions in `Externals.cfg`, and the `manage_externals/checkout_externals` script downloads them.
- **CAM** is the atmosphere component. **WACCM** is CAM with the model top raised to about 140 km (70 levels), which is needed here because alumina is released at 30-90 km.
- **CIME** is the set of scripts that create, configure, build and run a case (`create_newcase`, `case.setup`, `case.build`, `case.submit`).
- **CARMA** (Community Aerosol and Radiation Model for Atmospheres) is a *sectional* (bin) aerosol model that runs inside CAM. Instead of tracking a few bulk size modes, it splits particles into mass bins. Each bin is its own advected tracer in CAM. Bin `i` holds particles of mass `rmin_mass * vmrat^(i-1)`: with `vmrat = 2`, the particle volume doubles from one bin to the next, so the radius grows by `2^(1/3)` ≈ 1.26 per bin.
  - Inside CAM, CARMA handles coagulation, sedimentation (gravitational settling), dry deposition at the surface, wet deposition (rain-out), and vertical diffusion. CAM's dynamics handles large-scale transport.
  - CARMA comes as a generic engine (`src/physics/carma/base`, an external) plus *models* (`src/physics/carma/models/<name>`). A model says which particle types exist, how many bins, which processes are on, and how particles are emitted. You pick a model at build time with `-carma <name>` in `CAM_CONFIG_OPTS`.
  - At build time, `models/<name>/` comes before `base/` in the source search path, so a model can replace a base file by shipping a file with the same name. The alumina models use this (see `vertical.F90` below).

## 2. How the alumina model works

The model represents **involatile aluminium oxide (Al2O3, alumina) particles produced when spacecraft and rocket stages ablate during reentry**. It was built from CARMA's `meteor_smoke` model.

### 2.1 Particles and bins

| Property | Value |
|---|---|
| Groups / elements | 1 group, 1 element, both called `alumina`, short name `CRALUM` |
| Density | 3.95 g/cm³ |
| Shape | sphere |
| Composition | involatile (no growth, evaporation or nucleation) |
| Processes | coagulation (alumina-alumina), sedimentation, dry deposition, wet deposition (`solfac = 0.3`, `scavcoef = 0.1`), vertical diffusion, emission |
| Tracers in CAM | `CRALUM01`...`CRALUMnn` (mass mixing ratio, kg/kg, one per bin) |
| Radiation | radiatively active in RRTMG through per-bin optics files |

There are two kinds of alumina model:

- **`alumina`** (base model): 36 bins, `rmin = 0.5 nm`, `vmrat = 2`, so the bins span radii from about 0.5 nm to about 1.6 µm. It is meant for a realistic, altitude-dependent size distribution.
- **`alumina_d<D>um`** (particle-size sensitivity models): 16 bins, `vmrat = 2`, with `rmin` equal to the emission radius `D/2`. A monodisperse emission puts 100% of the mass in bin 1, and bins 2-16 (up to 32× the emission radius) catch particles that grow by coagulation. `D` is the diameter in µm; `p` stands for a decimal point. The available sizes are:

  `0p01, 0p1, 2, 6, 10, 14, 18, 25, 35, 45, 50, 100, 150, 200, 250`

  Each `alumina_d<D>um` directory has its own `carma_model_mod.F90` (only `NBIN`, `rmin` and the header comment differ from the base model). Its other files are symlinks to `../alumina/`, so a fix in the base model applies to all of them.

### 2.2 Emissions: discrete events

Emissions are not a continuous flux. They are **point releases** (Dirac pulses in time) read from a netCDF event file (`carma_emis_file`), handled in `alumina_events_mod.F90`:

- Each row of the file is one release: `event_id`, `date` (YYYYMMDD), `datesec` (seconds of the day), `lat`, `lon` (0-360), `altitude` (km above ground) and `mass` (kg of Al2O3). One reentry is normally several rows, one per altitude band.
- At startup, each row is assigned to the **nearest physics column** (great-circle distance), and its mass is split across the CARMA bins using the size distribution at its altitude (section 2.3).
- During the run, a row fires in the timestep where `start of step <= row time < end of step`. All of its mass goes into that one column and the model layer that contains its altitude, within that single timestep. Each row fires exactly once.
- `carma_emis_scale` multiplies all masses.
- With `carma_emis_cyclic = .true.`, rows are matched on month, day and time of day only, so the same year of events repeats every model year. The file must then hold one year of events. Rows never fire before the run start date.
- At startup the log reports the number of rows, total mass, the receiving column of each row, and how many rows are dated before the run start (those never fire in a non-cyclic run).

### 2.3 Particle size distribution (PSD)

`alumina_psd_mod.F90` decides how each row's mass is split among the bins:

- `carma_emis_psd_type = 'file'` (default): reads `carma_emis_psd_file`, which holds mass fractions on its own bins (`r_edge`, in nm) at a set of altitudes (`psd_alt`, km). At startup each altitude row is remapped onto the CARMA bins by overlap in ln(r), assuming mass is uniform in ln(r) inside each input bin. At emission, the PSD is linearly interpolated to the row's altitude (and held constant outside the range of `psd_alt`).
- `carma_emis_psd_type = 'lognormal'`: a single lognormal at all altitudes, with number median radius `carma_emis_rmode` (nm) and geometric standard deviation `carma_emis_sigma`.
- Mass that falls outside the CARMA bin range is "clipped". If more than 1% of the mass is clipped, the run stops, unless `carma_emis_psd_allow_clip = .true.` (then it only warns). After clipping, each row is renormalized to 1.
- For the size models, the PSD files are "monodisperse": one very narrow input bin at radius `D/2`, so all the mass lands in CARMA bin 1. Note that the particles then have the size of bin 1, not exactly `D/2`.

### 2.4 Sedimentation override (`vertical.F90`)

CARMA normally skips sedimentation in any column whose peak number concentration is below `FEW_PC = 1e-6 /cm³`. Large alumina particles can carry a lot of mass at number concentrations well below that, so they would never fall. `models/alumina/vertical.F90` is a copy of CARMA's `base/vertical.F90` with that test removed (search for `ALUMINA OVERRIDE`). Because model files come before `base/` in the build path, it replaces the base version for all alumina models. Coagulation and the other processes keep the `FEW_PC` test.

### 2.5 Radiation and the radiative effect

Each CARMA bin needs an optics file (`<model>_CRALUMnn_rrtmg.nc`) with its extinction, single-scattering albedo and asymmetry factor in each RRTMG band. These come from Mie theory using alumina refractive indices that are hard-coded in `CARMA_DefineModel` (16 LW and 14 SW band centres).

- **Base `alumina` model:** CARMA writes its 36 optics files into the run directory at startup (`carma_do_optics = .true.`).
- **`alumina_d<D>um` models:** CARMA's default Mie routine (`miess`) fails for alumina particles larger than about 175 µm. So these models use optics files **pre-built** with `components/cam/tools/carma_optics_offline` (patched Bohren & Huffman `bhmie`), stored in `$DIN_LOC_ROOT/atm/cam/physprops/alumina/`, with `carma_do_optics = .false.`.

CAM's radiation can be called several times per step with different aerosol lists. `build-namelist` sets them up like this:

| List | Contents | Purpose |
|---|---|---|
| `rad_climate` | standard aerosols + alumina bins **if** `carma_alumina_rad_feedback = .true.`; standard aerosols only if `.false.` | the radiation that actually heats and cools the model |
| `rad_diag_1` | standard aerosols + alumina bins | diagnostic call (output `*_d1`) |
| `rad_diag_2` | standard aerosols only | diagnostic call (output `*_d2`) |

So **the alumina radiative effect is `X_d1 - X_d2`** (for example `FSNT_d1 - FSNT_d2`), whether or not feedback is on. Don't difference `X_d1` against the climate call `X`: MAM recomputes aerosol water uptake for diagnostic calls, so those differ even with identical lists.

- `carma_alumina_rad_feedback = .true.`: alumina changes the simulated climate.
- `carma_alumina_rad_feedback = .false.` ("no feedback", `nofb`): alumina is transported and deposited but does not affect the climate, so the weather stays close to a run without alumina. The radiative effect is still diagnosed by `X_d1 - X_d2`.

## 3. What was changed in CAM to add it

All changes are in the CAM repository. The commits on `ashajain/carma_aluminum_model`, on top of NCAR's `cam_cesm2_1_rel` (commit `a03b84b7c`), are:

| Commit | Change |
|---|---|
| Register CARMA alumina model in the CAM build system | `bld/configure`, `bld/config_files/definition.xml` |
| Add CARMA alumina model for spacecraft-reentry Al2O3 | `src/physics/carma/models/alumina/*`, namelist entries |
| Make CARMA alumina radiatively active in RRTMG | `bld/build-namelist` (`rad_climate`, `rad_diag_1`, `rad_diag_2`), `carma_alumina_rad_feedback` |
| Add CARMA alumina particle size sensitivity models | `src/physics/carma/models/alumina_d<D>um/` |
| Use pre-built optics files for the alumina size models | `build-namelist` points at `$DIN_LOC_ROOT/atm/cam/physprops/alumina/` |
| Add offline CARMA optics tool for the alumina models | `tools/carma_optics_offline/` |
| Always sediment alumina particles | `models/alumina/vertical.F90` |
| Add `carma_emis_cyclic` option | repeat the emission events every model year |

File by file:

**New CARMA model code** (`src/physics/carma/models/alumina/`):
- `carma_model_mod.F90`: defines the group, element, bins, refractive indices and coagulation (`CARMA_DefineModel`), sets up the PSD and events (`CARMA_InitializeModel`), and adds the emission tendency (`CARMA_EmitParticle`).
- `carma_model_flags_mod.F90`: the alumina namelist variables (`carma_model_nl` group) and their MPI broadcast.
- `alumina_events_mod.F90`: reads the event file, maps rows to columns, fires them at the right timestep.
- `alumina_psd_mod.F90`: reads or builds the PSD and remaps it onto the CARMA bins.
- `vertical.F90`: the sedimentation override (section 2.4).
- `models/alumina_d<D>um/`: one directory per size; only `carma_model_mod.F90` is a real file, the rest are symlinks.

**Build and namelist system** (`bld/`):
- `configure`: accepts `-carma alumina` and `-carma alumina_d<D>um`; sets the number of advected tracers (36 or 16); raises `max_n_rad_cnst` to 64 because the radiation lists become longer than the default limit of 30.
- `config_files/definition.xml`: adds the model names to the valid values of `carma`.
- `build-namelist`: sets the default CARMA process flags for alumina, builds the `rad_climate` / `rad_diag_1` / `rad_diag_2` lists, and points the size models at the pre-built optics files.
- `namelist_files/namelist_definition.xml`: defines the new namelist variables (`carma_emis_scale`, `carma_emis_cyclic`, `carma_emis_psd_type`, `carma_emis_psd_file`, `carma_emis_psd_allow_clip`, `carma_emis_rmode`, `carma_emis_sigma`, `carma_alumina_rad_feedback`; `carma_emis_file` is reused from `meteor_smoke`).
- `namelist_files/namelist_defaults_cam.xml`: default event and PSD files per model. **These default paths (`atm/waccm/emis/alumina_*`) are not in NCAR's inputdata, so always set `carma_emis_file` and `carma_emis_psd_file` in `user_nl_cam`.**

**Tool** (`tools/carma_optics_offline/`): builds a stand-alone driver from CARMA's own source and writes the RRTMG optics files for the size models. See its `README.md`.

Nothing outside CAM was changed. The other components (CLM, CICE, CIME, etc.) are the standard CESM2.1.5 versions.

## 4. Installation: getting the code

These steps set up a fresh CESM2.1.5 checkout that uses the alumina CAM branch. Your machine must already be able to run CESM2.1 (a CIME machine entry with compilers, MPI, netCDF and a batch system). If it can't, set that up first (see the CIME porting guide), or ask someone who has already ported CESM to your cluster.

### 4.1 Clone CESM2.1.5

```bash
mkdir -p ~/cesm && cd ~/cesm
git clone -b release-cesm2.1.5 https://github.com/ESCOMP/CESM.git cesm2.1.5_alumina
cd cesm2.1.5_alumina
```

### 4.2 Point `Externals.cfg` at the alumina CAM branch

Edit the `[cam]` block of `Externals.cfg` so that it reads:

```ini
[cam]
branch = ashajain/carma_aluminum_model
protocol = git
repo_url = https://github.com/ashajain652/CAM_carma19
local_path = components/cam
externals = Externals_CAM.cfg
required = True
```

(Delete the `tag = ...` line that was there; a block must have either `tag` or `branch`, not both.) Leave the other blocks as they are.

### 4.3 Check out all the components

```bash
./manage_externals/checkout_externals
./manage_externals/checkout_externals --status --verbose   # every line should be clean (no 'e' or 'M')
```

This also checks out CAM's own externals (`Externals_CAM.cfg`), including the CARMA engine into `components/cam/src/physics/carma/base` (`carma3_49_rel`), which is fetched with **svn**. You need the `svn` command on your path. If the svn checkout fails, see [Troubleshooting](#10-troubleshooting).

Check that you got the alumina code:

```bash
cd components/cam
git log --oneline -10       # should list the alumina commits
ls src/physics/carma/models | grep alumina
cd ../..
```

To pick up later updates to the branch:

```bash
cd components/cam
git fetch origin
git checkout ashajain/carma_aluminum_model     # or: git pull, if you're on the branch
```

### 4.4 Get the emission and PSD files

```bash
cd ~/cesm
git clone https://github.com/ashajain652/alumina_emissions.git
```

This holds the event files, the PSD files and the scripts that make them (section 5). Ask Asha for access if the repository is private.

### 4.5 Get or build the optics files (size models only)

The `alumina_d<D>um` models read their optics from `$DIN_LOC_ROOT/atm/cam/physprops/alumina/alumina_d<D>um_CRALUMnn_rrtmg.nc` (16 files per size). `DIN_LOC_ROOT` is your machine's CESM inputdata folder (`./xmlquery DIN_LOC_ROOT` in any case).

Either copy them from Asha, or build them yourself (needs `gfortran`, and Python 3 with `numpy` and `netCDF4`):

```bash
cd ~/cesm/cesm2.1.5_alumina/components/cam/tools/carma_optics_offline
./build.sh
# all sizes:
./make_alumina_optics.py --all-sizes --outdir $DIN_LOC_ROOT/atm/cam/physprops/alumina
# or one size:
./make_alumina_optics.py alumina_d50um --outdir $DIN_LOC_ROOT/atm/cam/physprops/alumina
```

**Always pass `--outdir`**; the script's built-in default is a path on Asha's laptop.

The base `alumina` model doesn't need this step: it makes its optics at runtime.

## 5. Input files

| File | Namelist variable | Where it comes from |
|---|---|---|
| Event file | `carma_emis_file` | `alumina_emissions/alumina_2025_reentry_events.nc` (1223 reentries in 2025, 192,194.64 kg Al2O3 in total, 30-90 km in 5 km bands, dates shifted to model year 0001). A small test file is `alumina_test_events.nc`. |
| PSD file, base model | `carma_emis_psd_file` | `alumina_emissions/alumina_test_psd.nc` |
| PSD file, size models | `carma_emis_psd_file` | `alumina_emissions/psd_sensitivitytest/alumina_mono_d<D>um_psd.nc` |
| Optics, size models | (set by build-namelist) | `$DIN_LOC_ROOT/atm/cam/physprops/alumina/` (section 4.5) |

**Event file format** (netCDF, dimension `row`):

```
event_id(row)  int    event identifier (shared by the rows of one event)
date(row)      int    YYYYMMDD in the model calendar (noleap)
datesec(row)   int    seconds of the day, 0-86399
lat(row)       real   degrees north
lon(row)       real   degrees east, 0-360
altitude(row)  real   km above ground
mass(row)      real   kg of Al2O3
```

**The event dates must fall inside your run.** The 2025 file is dated year 0001, so runs should start at `0001-01-01`; with `carma_emis_cyclic = .true.` the events then repeat every model year. Rows dated before the start date never fire (the log reports how many).

To make your own inputs, use the scripts in `alumina_emissions` (each has usage notes in its header; they need `numpy` and `netCDF4`):

```bash
# reentry inventory -> emissions CSV (one row per reentry and 5 km band)
./reentry_csv_to_alumina_emissions.py --reentries my_reentries.csv -o my_emissions.csv
# emissions CSV (id,start_date,lat,long,altitude_km,emitted_mass_kg) -> event file
#   --year-offset shifts dates into the model calendar, e.g. -2024 maps 2025 to year 0001
./csv_to_alumina_events.py my_emissions.csv my_events.nc --year-offset -2024
./check_alumina_events.py my_events.nc            # prints rows and total mass

# PSD CSV (altitude_km,bin,d_low_um,d_center_um,d_high_um,dNdlnr_cm-3) -> PSD file
./csv_to_alumina_psd.py my_psd.csv my_psd.nc
# monodisperse PSD for a size model: diameter (um), and the model's rmin (cm) and NBIN
./make_monodisperse_psd.py 50 --rmin 2.5e-3 --nbin 16 -o alumina_mono_d50um_psd.nc
```

## 6. Running a case, step by step

The example below is the setup used for the particle-size experiments: WACCM with specified chemistry and a 2010 climatology (`FWsc2010climo`), at 1.9°×2.5° resolution, with 50 µm particles. Replace `<machine>` with your CIME machine name and the paths with your own.

### Step 1: create the case

```bash
cd ~/cesm/cesm2.1.5_alumina/cime/scripts
./create_newcase --case ~/cesm/cases/FWsc2010climo.alum50um \
                 --compset FWsc2010climo --res f19_f19_mg17 \
                 --machine <machine> --run-unsupported
cd ~/cesm/cases/FWsc2010climo.alum50um
```

`FWsc2010climo` = `2010_CAM60%WCSC_CLM50%SP_CICE%PRES_DOCN%DOM_MOSART_CISM2%NOEVOLVE_SWAV`: CAM6 physics, WACCM with specified chemistry (MAM4 aerosols), prescribed sea surface temperature and sea ice, and 2010 forcing repeated every year.

### Step 2: switch on CARMA

Add the CARMA model to the CAM build options. **This must be done before `case.setup`/`case.build`; changing it later needs a clean rebuild.**

```bash
./xmlchange --append CAM_CONFIG_OPTS="-carma alumina_d50um"
./xmlquery CAM_CONFIG_OPTS
#  -> -phys cam6 -age_of_air_trcs -chem waccm_sc_mam4 -carma alumina_d50um
```

Use `-carma alumina` for the 36-bin base model, or any `alumina_d<D>um` from section 2.1.

### Step 3: set the run length, start date and processors

```bash
./xmlchange RUN_TYPE=startup,RUN_STARTDATE=0001-01-01   # matches the year-0001 event dates
./xmlchange STOP_OPTION=nmonths,STOP_N=4                # one job = 4 model months
./xmlchange REST_OPTION=nmonths,REST_N=4                # restart files at the end of each job
./xmlchange RESUBMIT=2                                  # 3 jobs in total = 1 year
./xmlchange JOB_WALLCLOCK_TIME=12:00:00
./xmlchange NTASKS=48                                   # all components; adjust to your nodes
```

For scale: at f19 with 48 MPI tasks, this model runs at about 140-200 s per model day, so about 7 h per 4-month job. Pick `STOP_N` so that one job fits in your queue's wall-clock limit.

### Step 4: set up the case

```bash
./case.setup
```

### Step 5: edit `user_nl_cam`

A minimal `user_nl_cam` for a size model (see section 7 for every option and section 8 for output):

```fortran
! Alumina emissions and size distribution (absolute paths)
carma_emis_file     = '/home/<you>/cesm/alumina_emissions/alumina_2025_reentry_events.nc'
carma_emis_psd_file = '/home/<you>/cesm/alumina_emissions/psd_sensitivitytest/alumina_mono_d50um_psd.nc'
! Repeat the year-0001 events every model year
carma_emis_cyclic   = .true.
! .true.: alumina affects the climate; .false.: diagnose the radiative effect only
carma_alumina_rad_feedback = .true.

! Output: monthly means
nhtfrq = 0, 0, 0
mfilt  = 1, 1, 1
ndens  = 2, 2, 1            ! h2 (radiation) in double precision: alumina flux differences are tiny
fincl2 = 'CRALUMMR', 'CRALUMND', 'CRALUMRE', 'CRALUMVM',
         'CRALUM01EM', 'CRALUM01DD', 'CRALUM01SW', 'CRALUM01WD'      ! repeat for bins 02-16
fincl3 = 'FSNT', 'FSNT_d1', 'FSNT_d2', 'FLNT', 'FLNT_d1', 'FLNT_d2',
         'FSNS', 'FSNS_d1', 'FSNS_d2', 'FLNS', 'FLNS_d1', 'FLNS_d2',
         'QRS', 'QRS_d1', 'QRS_d2', 'QRL', 'QRL_d1', 'QRL_d2'
```

The per-bin tracers `CRALUMnn` are not written by default: add them to a `fincl` list (as `'CRALUM01:I'` etc. for instantaneous values), or set `history_carma = .true.` to put them and the group fields on h0. A complete, working `user_nl_cam` from the size experiments is in [Appendix A](#appendix-a-full-user_nl_cam-from-the-size-experiments).

Optional: to start from a spun-up atmosphere instead of NCAR's default initial file, add `ncdata = '/path/to/cam.i.nc'` (and `finidat` in `user_nl_clm`).

### Step 6: check the namelists

```bash
./preview_namelists
grep -i carma CaseDocs/atm_in | head -40      # carma_model, carma_emis_*, carma_do_* flags
grep -A20 rad_climate CaseDocs/atm_in         # CRALUM bins present (feedback on) or absent (off)
grep -A20 rad_diag_1  CaseDocs/atm_in         # always lists the CRALUM bins
./check_input_data                             # all input files found
```

`preview_namelists` stops with an error if a `user_nl_cam` variable is misspelled or a file doesn't exist.

### Step 7: build

```bash
./case.build
```

A successful build ends with `MODEL BUILD HAS FINISHED SUCCESSFULLY`. Changes to `user_nl_cam` do **not** need a rebuild; changes to `CAM_CONFIG_OPTS` or the Fortran source do (`./case.build --clean-all` then `./case.build`).

### Step 8: submit

```bash
./case.submit
```

With `RESUBMIT=2`, each job resubmits the next one when it finishes. Follow progress with:

```bash
cat CaseStatus
tail -f $(./xmlquery --value RUNDIR)/atm.log.*      # 'nstep' lines show model time
```

### Step 9: check the first job

In the run directory (`./xmlquery RUNDIR`):

```bash
grep -i "alumina_events_init" atm.log.*    # rows, total mass, carma_emis_scale, carma_emis_cyclic
grep -i "clipped"             atm.log.*    # PSD mass clipped outside the bins (should be ~0)
grep -i "WARNING"             atm.log.*    # e.g. rows dated before the run start
```

When finished, output is in `$DOUT_S_ROOT/atm/hist/` (`./xmlquery DOUT_S_ROOT`) as `<case>.cam.h0.*`, `h1`, ... files.

### Continuing or changing a run

- To extend a run: `./xmlchange CONTINUE_RUN=TRUE,RESUBMIT=<n>` and `./case.submit`.
- **The output frequency (`nhtfrq`, `mfilt`) is read back from the restart file**, so changing it in `user_nl_cam` has no effect on a continued run. Start a new run (startup) or a branch run from a restart instead.
- To make a copy of a case with a different setting, use `create_clone`:

  ```bash
  cd ~/cesm/cesm2.1.5_alumina/cime/scripts
  # same CARMA model, different namelist only: reuse the executable
  ./create_clone --case ~/cesm/cases/FWsc2010climo.alum50um.nofb \
                 --clone ~/cesm/cases/FWsc2010climo.alum50um --keepexe
  # different CARMA model: build a new executable (no --keepexe)
  ./create_clone --case ~/cesm/cases/FWsc2010climo.alum10um \
                 --clone ~/cesm/cases/FWsc2010climo.alum50um
  ```

  With `--keepexe`, don't rebuild or clean the original case while the clone uses its executable. Without it, check that `EXEROOT` points at the new case's own `bld` folder before building.

### Suggested set of runs

| Run | How |
|---|---|
| Alumina, with feedback | as above, `carma_alumina_rad_feedback = .true.` |
| Alumina, no feedback | clone `--keepexe`, set `carma_alumina_rad_feedback = .false.`, and remove `ODV_CRALUMnn` from `fincl` (section 8) |
| Reference without alumina | same compset and dates, no `-carma` |
| Noise estimate | clones with `pertlim = 1.0e-14`, `2.0e-14` in `user_nl_cam`, to measure run-to-run (weather) variability |

## 7. `user_nl_cam` reference

Alumina-specific variables (namelist group `carma_model_nl`):

| Variable | Default | Meaning |
|---|---|---|
| `carma_emis_file` | test file (not in NCAR inputdata) | **Set this.** Event file, absolute path. |
| `carma_emis_psd_file` | per-model file (not in NCAR inputdata) | **Set this** when `carma_emis_psd_type = 'file'`. Absolute path. |
| `carma_emis_scale` | `1.0` | Multiplier on all event masses. |
| `carma_emis_cyclic` | `.false.` | `.true.` repeats the events every model year (file must hold one year). |
| `carma_emis_psd_type` | `'file'` | `'file'` or `'lognormal'`. |
| `carma_emis_psd_allow_clip` | `.false.` | `.true.` only warns if >1% of PSD mass falls outside the bins. |
| `carma_emis_rmode` | `2.0` | Lognormal number median radius (nm), `'lognormal'` only. |
| `carma_emis_sigma` | `1.5` | Lognormal geometric standard deviation, `'lognormal'` only. |
| `carma_alumina_rad_feedback` | `.true.` | Include alumina in `rad_climate`. |

General CARMA flags (group `carma_nl`), set by build-namelist for alumina; you normally leave them alone:

| Variable | Alumina default | Meaning |
|---|---|---|
| `carma_do_emission` | `.true.` | emissions on |
| `carma_do_coag` | `.true.` | coagulation |
| `carma_do_vtran` | `.true.` | sedimentation |
| `carma_do_drydep` | `.true.` | dry deposition |
| `carma_do_wetdep` | `.true.` | wet deposition |
| `carma_do_vdiff` | `.true.` | vertical diffusion |
| `carma_do_fixedinit` | `.true.` | initialize from fixed reference state |
| `carma_do_optics` | `.true.` (base) / `.false.` (size models) | write optics files at runtime |

Setting `rad_climate`, `rad_diag_1` or `rad_diag_2` yourself overrides what `carma_alumina_rad_feedback` sets up, so avoid that unless you know what you need.

Useful general CAM settings:

| Variable | Example | Meaning |
|---|---|---|
| `nhtfrq` | `0` monthly, `-24` daily, `-72` every 3 days, `1` every step | output frequency per history file (h0, h1, ...) |
| `mfilt` | `1`, `30`, ... | time samples per file |
| `ndens` | `2` float, `1` double | precision per file |
| `fincl1`...`fincl6` | field lists | fields added to each history file; `:I` = instantaneous instead of mean |
| `empty_htapes` | `.true.` | drop the default h0 fields so only your `fincl1` fields are written |
| `ncdata` | path | atmosphere initial condition |
| `pertlim` | `1.0e-14` | tiny random temperature perturbation at startup (ensembles) |

## 8. Output fields and analysis tips

`nn` is the bin number (01-36 for `alumina`, 01-16 for the size models).

| Field | Units | Meaning |
|---|---|---|
| `CRALUMnn` | kg/kg | mass mixing ratio of bin `nn` |
| `CRALUMMR` | kg/kg | total alumina mass mixing ratio |
| `CRALUMND` | #/cm³ | number density |
| `CRALUMMD` | g/cm³ | mass density |
| `CRALUMRE` | µm | effective radius |
| `CRALUMVM` | m/s | fall velocity |
| `CRALUMnnEM` | kg/kg/s | emission tendency (3-D) |
| `CRALUMnnTC` | kg/kg/s | total CARMA tendency (3-D) |
| `CRALUMnnDD` | kg/m²/s | sedimentation + dry deposition flux at the surface |
| `CRALUMnnSW` | kg/m²/s | wet deposition flux at the surface |
| `CRALUMnnWD` | kg/kg/s | wet deposition tendency (3-D) |
| `FSNT_d1`, `FSNT_d2`, `FLNT_d1`, ... | W/m² | radiation diagnostics with (`d1`) and without (`d2`) alumina |
| `QRS_d1`, `QRL_d1`, ... | K/s | SW/LW heating rates with/without alumina |
| `ODV_CRALUMnn` | 1 | visible optical depth of bin `nn` |

Tips:

- **Radiative effect** = `X_d1 - X_d2`. These differences are tiny (µW/m² for this inventory), so write them in double precision (`ndens = 1` for that file).
- **`ODV_CRALUMnn` only exists when alumina is in `rad_climate`.** With `carma_alumina_rad_feedback = .false.`, asking for it in `fincl` makes the run stop at startup.
- **Burden:** integrate `CRALUMnn` over the column using the layer air mass (output `PS` with the hybrid coefficients, or the `MASS` field). Instantaneous (`:I`) values at the end of each output period make a mass budget easier to close.
- **Large particles fall out in days.** For particles larger than about 10 µm, monthly means hide most of the signal; use daily (`nhtfrq = -24`) or 3-day output, and keep the field lists short, as 3-D daily output is large.
- **Noise:** the month-to-month change in net TOA flux from internal variability is several W/m², far larger than the alumina effect. Only `X_d1 - X_d2` within one run is clean; perturbed-minus-reference differences between runs are dominated by weather noise unless you average an ensemble.

## 9. Adding a new particle-size model

Example: a 75 µm diameter model, `alumina_d75um` (`rmin` = 37.5 µm = 3.75e-3 cm).

1. Copy an existing size model, keeping the symlinks:
   ```bash
   cd components/cam/src/physics/carma/models
   cp -a alumina_d50um alumina_d75um
   ```
2. In `alumina_d75um/carma_model_mod.F90`, set `rmin = 3.7500e-03_f` and update the header comment.
3. In `components/cam/bld/config_files/definition.xml`, add `alumina_d75um` to the `valid_values` of `carma`. (`bld/configure` already accepts any name matching `alumina_d[0-9p]+um` and sets 16 bins.)
4. In `bld/namelist_files/namelist_defaults_cam.xml`, add `carma_emis_file` and `carma_emis_psd_file` lines for `carma="alumina_d75um"`, following the existing ones.
5. Add `75` to `SIZES` in `tools/carma_optics_offline/make_alumina_optics.py`, then make its optics:
   ```bash
   ./make_alumina_optics.py alumina_d75um --outdir $DIN_LOC_ROOT/atm/cam/physprops/alumina
   ```
6. Make its PSD file:
   ```bash
   ./make_monodisperse_psd.py 75 --rmin 3.75e-3 --nbin 16 -o psd_sensitivitytest/alumina_mono_d75um_psd.nc
   ```
7. Create a case with `-carma alumina_d75um` and check in `atm.log` that 100% of the mass goes into bin 1 and nothing is clipped.

If you change the number of bins, also update `bld/configure` (`$carma_nadv`) and `bld/build-namelist` (`$nbin`), which assume 16 for the size models.

## 10. Troubleshooting

| Problem | Fix |
|---|---|
| `checkout_externals` fails on `carma` (svn) | Install `svn`, or accept the server certificate by running the svn checkout by hand once. If the svn server can't be reached, copy `components/cam/src/physics/carma/base` from a working checkout (Asha's) into yours. |
| `preview_namelists` / `check_input_data` can't find `atm/waccm/emis/alumina_*.nc` | Set `carma_emis_file` and `carma_emis_psd_file` to absolute paths in `user_nl_cam`. |
| Run stops: optics file `alumina_d<D>um_CRALUMnn_rrtmg.nc` not found | Build or copy the optics files into `$DIN_LOC_ROOT/atm/cam/physprops/alumina/` (section 4.5). |
| Run stops at startup on an `ODV_CRALUMnn` field | Remove `ODV_*` from `fincl` when `carma_alumina_rad_feedback = .false.`. |
| Run stops: PSD mass clipped | The PSD extends beyond the model's bins. Use a model whose bins cover it, fix the PSD, or set `carma_emis_psd_allow_clip = .true.` if a little loss is acceptable. |
| No alumina in the output | Check the event dates against the run dates (`grep alumina_events_init atm.log.*`). Non-cyclic rows dated before `RUN_STARTDATE` never fire. |
| `carma_emis_cyclic needs all rows in one year` | With cyclic on, the event file must hold a single year. |
| Changing `nhtfrq` has no effect | Output frequency comes from the restart file on a continued run; start a new or branch run. |
| Changed `-carma` but nothing changed | `./case.build --clean-all` then `./case.build`. |
| `Undefined env var 'LD_LIBRARY_PATH'` from CIME scripts (seen on MIT SuperCloud) | `export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}` before running them. |

---

## Appendix A: full `user_nl_cam` from the size experiments

This is the `user_nl_cam` of the 50 µm no-feedback case. For a feedback case, set `carma_alumina_rad_feedback = .true.` and you may add `'ODV_CRALUM01'`...`'ODV_CRALUM16'` to `fincl4`.

```fortran
! CARMA alumina inputs: 2025 reentries, dates shifted to model year 1. No emission scaling.
carma_emis_file     = '/path/to/alumina_emissions/alumina_2025_reentry_events.nc'
carma_emis_psd_file = '/path/to/alumina_emissions/psd_sensitivitytest/alumina_mono_d50um_psd.nc'
! Repeat the year-0001 events every model year (matched on MMDD + seconds).
carma_emis_cyclic   = .true.

! No radiative feedback: alumina is left out of rad_climate. It is still in rad_diag_1
! (not rad_diag_2), so X_d1 - X_d2 is its radiative effect.
carma_alumina_rad_feedback = .false.

empty_htapes = .true.

! h0: monthly means of the basic state.
! h1: monthly means of per-bin emission (EM) and losses (DD, SW, WD).
! h2: end-of-month instantaneous CRALUM mixing ratios, for the burden.
! h3: monthly means of the alumina radiative effect (double precision).
nhtfrq = 0, 0, 0, 0
mfilt  = 1, 1, 1, 1
ndens  = 2, 2, 2, 1
fincl1 = 'PS', 'T', 'U', 'V', 'Q', 'Z3', 'OMEGA',
        'TS', 'TREFHT', 'PRECT', 'CLDTOT', 'TROP_P', 'PMID'
fincl2 = 'PS', 'MASS', 'AREA', 'CRALUMVM', 'CRALUM01EM', 'CRALUM02EM', 'CRALUM03EM', 'CRALUM04EM', 'CRALUM05EM',
        'CRALUM06EM', 'CRALUM07EM', 'CRALUM08EM', 'CRALUM09EM', 'CRALUM10EM', 'CRALUM11EM',
        'CRALUM12EM', 'CRALUM13EM', 'CRALUM14EM', 'CRALUM15EM', 'CRALUM16EM', 'CRALUM01DD',
        'CRALUM02DD', 'CRALUM03DD', 'CRALUM04DD', 'CRALUM05DD', 'CRALUM06DD', 'CRALUM07DD',
        'CRALUM08DD', 'CRALUM09DD', 'CRALUM10DD', 'CRALUM11DD', 'CRALUM12DD', 'CRALUM13DD',
        'CRALUM14DD', 'CRALUM15DD', 'CRALUM16DD', 'CRALUM01SW', 'CRALUM02SW', 'CRALUM03SW',
        'CRALUM04SW', 'CRALUM05SW', 'CRALUM06SW', 'CRALUM07SW', 'CRALUM08SW', 'CRALUM09SW',
        'CRALUM10SW', 'CRALUM11SW', 'CRALUM12SW', 'CRALUM13SW', 'CRALUM14SW', 'CRALUM15SW',
        'CRALUM16SW', 'CRALUM01WD', 'CRALUM02WD', 'CRALUM03WD', 'CRALUM04WD', 'CRALUM05WD',
        'CRALUM06WD', 'CRALUM07WD', 'CRALUM08WD', 'CRALUM09WD', 'CRALUM10WD', 'CRALUM11WD',
        'CRALUM12WD', 'CRALUM13WD', 'CRALUM14WD', 'CRALUM15WD', 'CRALUM16WD'
fincl3 = 'CRALUM01:I', 'CRALUM02:I', 'CRALUM03:I', 'CRALUM04:I', 'CRALUM05:I', 'CRALUM06:I',
        'CRALUM07:I', 'CRALUM08:I', 'CRALUM09:I', 'CRALUM10:I', 'CRALUM11:I', 'CRALUM12:I',
        'CRALUM13:I', 'CRALUM14:I', 'CRALUM15:I', 'CRALUM16:I', 'PS:I', 'MASS:I'
fincl4 = 'FSNT', 'FSNT_d1', 'FSNT_d2', 'FLNT', 'FLNT_d1', 'FLNT_d2',
        'FSNS', 'FSNS_d1', 'FSNS_d2', 'FLNS', 'FLNS_d1', 'FLNS_d2',
        'FSUTOA', 'FSUTOA_d1', 'FSUTOA_d2', 'FLUT', 'FLUT_d1', 'FLUT_d2',
        'FSNR', 'FSNR_d1', 'FSNR_d2', 'FLNR', 'FLNR_d1', 'FLNR_d2',
        'FSNTC', 'FSNTC_d1', 'FSNTC_d2', 'FLNTC', 'FLNTC_d1', 'FLNTC_d2',
        'QRS', 'QRS_d1', 'QRS_d2', 'QRL', 'QRL_d1', 'QRL_d2',
        'AEROD_v', 'AODvstrt'

! Optional: start from a spun-up state (CESM2.1 FWsc2010climo year 3, regridded f09 -> f19)
! ncdata = '/path/to/f.e21.FWsc2010climo.f09_f09_mg17.cesm2.1-exp011.001_v2.cam.i.0003-01-01-00000.regrid_1.9x2.5.nc'
```

For the 36-bin base `alumina` model, extend the per-bin lists to `CRALUM36`.
