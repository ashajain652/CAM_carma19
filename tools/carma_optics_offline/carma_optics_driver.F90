!! Offline driver for CARMA's own Mie code: computes the RRTMG optical properties
!! that CARMA_CreateOpticsFile (cam/carma_intr.F90) writes for each bin of a
!! one-group, one-element CARMA model, without running CAM.
!!
!! It mirrors carma_register + CARMA_CreateOpticsFile:
!!   - band centres: get_{lw,sw}_spectral_boundaries('cm'), dwave = hi - lo,
!!     wave = lo + dwave/2
!!   - CARMA_Create, CARMAGROUP_Create (same arguments as the alumina models,
!!     plus the chosen Mie routine), CARMAELEMENT_Create, CARMA_Initialize
!!   - per bin: getwetr at mie_rh(1), then mie() at every band centre and the
!!     same unit conversions (per unit dry mass, mks)
!!
!! Input: namelist file given as the first command line argument (&optics_nl).
!! Output (stdout): lines tagged BIN, LW, SW with full precision values; parsed
!! by make_alumina_optics.py, which writes the netCDF files.
program carma_optics_driver

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carma_mod
  use radconstants, only : nswbands, nlwbands, get_lw_spectral_boundaries_cm, get_sw_spectral_boundaries_cm
  use wetr, only         : getwetr

  implicit none

  integer, parameter :: MAXBIN = 64
  integer, parameter :: I_GRP = 1, I_ELEM = 1, I_COMP = 1

  ! namelist
  integer            :: nbin      = 0
  real(kind=f)       :: rmin      = 0._f        ! cm
  real(kind=f)       :: vmrat     = 0._f
  real(kind=f)       :: rho_elem  = 0._f        ! g/cm3
  real(kind=f)       :: n_re(NWAVE) = 0._f
  real(kind=f)       :: n_im(NWAVE) = 0._f
  character(len=16)  :: mie_routine = 'bohren'  ! 'bohren' (bhmie) or 'toon' (miess)
  real(kind=f)       :: rh        = 0._f        ! mie_rh(1) of the models
  namelist /optics_nl/ nbin, rmin, vmrat, rho_elem, n_re, n_im, mie_routine, rh

  type(carma_type), target :: carma
  character(len=256) :: nlfile
  integer            :: rc, ibin, iwave, imiertn, unitn
  real(kind=f)       :: wave(NWAVE), dwave(NWAVE)
  logical            :: do_wave_emit(NWAVE)
  complex(kind=f)    :: refidx(NWAVE)
  real(kind=f), allocatable :: r(:), rlow(:), rup(:), rmass(:), rho(:)
  real(kind=f)       :: rwet, rhopwet, Qext, Qsca, asym

  call get_command_argument(1, nlfile)
  open(newunit=unitn, file=trim(nlfile), status='old')
  read(unitn, optics_nl)
  close(unitn)
  if (nbin < 1 .or. nbin > MAXBIN) stop 'bad nbin'

  select case (trim(mie_routine))
  case ('bohren')
    imiertn = I_MIERTN_BOHREN1983
  case ('toon')
    imiertn = I_MIERTN_TOON1981
  case default
    stop 'mie_routine must be bohren or toon'
  end select

  ! Band centres exactly as carma_register.
  call get_lw_spectral_boundaries_cm(wave(:nlwbands), dwave(:nlwbands))
  call get_sw_spectral_boundaries_cm(wave(nlwbands+1:), dwave(nlwbands+1:))
  dwave = dwave - wave
  wave  = wave + (dwave / 2._f)
  do_wave_emit(:nlwbands)  = .true.
  do_wave_emit(nlwbands+1:) = .true.
  do_wave_emit(nlwbands+1) = .false.
  do_wave_emit(NWAVE)      = .false.

  refidx(:) = cmplx(n_re(:), n_im(:), kind=f)

  call CARMA_Create(carma, nbin, 1, 1, 0, 0, NWAVE, rc, LUNOPRT=0, &
                    wave=wave, dwave=dwave, do_wave_emit=do_wave_emit)
  if (rc < 0) stop 'CARMA_Create failed'

  call CARMAGROUP_Create(carma, I_GRP, "alumina", rmin, vmrat, I_SPHERE, 1._f, .false., &
                         rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                         scavcoef=0.1_f, shortname="CRALUM", refidx=refidx, do_mie=.true., &
                         imiertn=imiertn)
  if (rc < 0) stop 'CARMAGROUP_Create failed'

  call CARMAELEMENT_Create(carma, I_ELEM, I_GRP, "alumina", rho_elem, I_INVOLATILE, I_COMP, rc, &
                           shortname="CRALUM")
  if (rc < 0) stop 'CARMAELEMENT_Create failed'

  ! Sets up the bins (setupbins), as in carma_register.
  call CARMA_Initialize(carma, rc, do_cnst_rlh=.true., do_print_init=.false.)
  if (rc < 0) stop 'CARMA_Initialize failed'

  allocate(r(nbin), rlow(nbin), rup(nbin), rmass(nbin), rho(nbin))
  call CARMAGROUP_Get(carma, I_GRP, rc, r=r, rlow=rlow, rup=rup, rmass=rmass)
  if (rc < 0) stop 'CARMAGROUP_Get failed'
  call CARMAELEMENT_Get(carma, I_ELEM, rc, rho=rho)
  if (rc < 0) stop 'CARMAELEMENT_Get failed'

  do iwave = 1, NWAVE
    write(*, '(a, i4, es25.16e3)') 'WAVE ', iwave, wave(iwave)
  end do

  do ibin = 1, nbin
    ! cgs values; converted to mks as CARMA_CreateOpticsFile does
    write(*, '(a, i4, 5es25.16e3)') 'BIN ', ibin, r(ibin), rlow(ibin), rup(ibin), rmass(ibin), rho(ibin)

    call getwetr(carma, I_GRP, rh, r(ibin), rwet, rho(ibin), rhopwet, rc)
    if (rc < 0) stop 'getwetr failed'

    do iwave = 1, NWAVE
      call mie(carma, carma%f_group(I_GRP)%f_imiertn, rwet, carma%f_wave(iwave), &
               carma%f_group(I_GRP)%f_nmon(ibin), carma%f_group(I_GRP)%f_df(ibin), &
               carma%f_group(I_GRP)%f_rmon, carma%f_group(I_GRP)%f_falpha, &
               carma%f_group(I_GRP)%f_refidx(iwave), Qext, Qsca, asym, rc)
      if (rc < 0) then
        write(*, '(a, 2i4)') 'MIEFAIL ', ibin, iwave
        stop 1
      end if

      if (iwave <= nlwbands) then
        write(*, '(a, 2i4, es25.16e3)') 'LW ', ibin, iwave, &
          (Qext - Qsca) * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f)
      else
        write(*, '(a, 2i4, 3es25.16e3)') 'SW ', ibin, iwave - nlwbands, &
          Qext * PI * (rwet * 1e-2_f)**2 / (rmass(ibin) * 1e-3_f), Qsca / Qext, asym
      end if
    end do
  end do

end program carma_optics_driver
