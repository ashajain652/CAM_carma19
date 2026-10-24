!! This module handles reading the namelist and provides access to some other flags
!! that control a specific CARMA model's behavior.
!!
!! By default the specific CARMA model does not have any unique namelist values. If
!! a CARMA model wishes to have its own namelist, then this file needs to be copied
!! from physics/cam to physics/model/<model_name> and the code needed to read in the
!! namelist values added there. This file will take the place of the one in
!! physics/cam.
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! This version is for the arsenic_msp model (reentry arsenic carried on meteoric
!! smoke). It is based upon the meteor_smoke version, with the emission scaling
!! (escale) options removed and the arsenic (carma_as_*) options added.
!!
!! @author  Chuck Bardeen
!! @version Mar-2011
!! @version Oct-2026 (arsenic_msp)
module carma_model_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics
  public carma_model_readnl                   ! read the carma model namelist


  ! Namelist flags
  !
  ! Create a public definition of any new namelist variables that you wish to have,
  ! and default them to an inital value.
  real(r8), public               :: carma_emis_total    = 16.0_r8  ! Total smoke mass emitted (kt/year)
  character(len=256), public     :: carma_emis_file     = 'atm/waccm/emis/meteor_smoke_kalashnikova.nc'   ! name of the smoke emission file
  real(r8), public               :: carma_as_emis_total = 4.0_r8   ! Global arsenic source, uniform in time and per unit area (Gg As/year)
  real(r8), public               :: carma_as_emis_alt   = 85.0_r8  ! Arsenic emission altitude (km above sea level)
  real(r8), public               :: carma_as_gamma      = 1.0_r8   ! Arsenic uptake (sticking) coefficient on smoke [0-1]
  logical, public                :: carma_as_do_uptake  = .true.   ! If .true. then do uptake of arsenic gas onto smoke
  logical, public                :: carma_as_do_fuchs   = .true.   ! If .true. then apply the Fuchs-Sutugin correction to the uptake
  logical, public                :: carma_as_do_binshift = .true.  ! If .true. then the arsenic taken up adds to the particle mass (bin shift)
  logical, public                :: carma_as_top_noflux = .true.   ! If .true. then no particle flux through the model top (DD is surface only)

contains


  !! Read the CARMA model runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Mar-2011
  subroutine carma_model_readnl(nlfile)

    ! Read carma namelist group.

    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    ! args

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars

    integer :: unitn, ierr

    ! read namelist for CARMA
    namelist /carma_model_nl/ &
      carma_emis_total, &
      carma_emis_file, &
      carma_as_emis_total, &
      carma_as_emis_alt, &
      carma_as_gamma, &
      carma_as_do_uptake, &
      carma_as_do_fuchs, &
      carma_as_do_binshift, &
      carma_as_top_noflux

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'carma_model_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, carma_model_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('carma_model_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(carma_emis_total,    1,                      mpir8,   0, mpicom)
    call mpibcast(carma_emis_file,     len(carma_emis_file),   mpichar, 0, mpicom)
    call mpibcast(carma_as_emis_total, 1,                      mpir8,   0, mpicom)
    call mpibcast(carma_as_emis_alt,   1,                      mpir8,   0, mpicom)
    call mpibcast(carma_as_gamma,      1,                      mpir8,   0, mpicom)
    call mpibcast(carma_as_do_uptake,  1,                      mpilog,  0, mpicom)
    call mpibcast(carma_as_do_fuchs,   1,                      mpilog,  0, mpicom)
    call mpibcast(carma_as_do_binshift, 1,                     mpilog,  0, mpicom)
    call mpibcast(carma_as_top_noflux, 1,                      mpilog,  0, mpicom)
#endif

  end subroutine carma_model_readnl

end module carma_model_flags_mod
