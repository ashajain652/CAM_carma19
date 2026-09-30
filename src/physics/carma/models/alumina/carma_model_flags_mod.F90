!! This module handles reading the namelist and provides access to some other flags
!! that control a specific CARMA model's behavior.
!!
!! This is the namelist for the alumina model, which represents involatile Al2O3
!! particles injected by spacecraft reentry ablation. Emissions are discrete point
!! events (dirac pulses in time) read from a netCDF event file, and the emitted mass
!! is split across the CARMA bins using an altitude dependent particle size
!! distribution (PSD) read from a netCDF file, or a lognormal for testing.
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen (template), alumina model
!! @version Sep-2026
module carma_model_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics
  public carma_model_readnl                   ! read the carma model namelist


  ! Namelist flags
  !
  ! Create a public definition of any new namelist variables that you wish to have,
  ! and default them to an inital value.
  character(len=256), public     :: carma_emis_file           = 'alumina_events.nc' ! event file (netCDF)
  real(r8), public               :: carma_emis_scale          = 1.0_r8              ! multiplier on all event masses
  character(len=16), public      :: carma_emis_psd_type       = 'file'              ! 'file' or 'lognormal'
  character(len=256), public     :: carma_emis_psd_file       = 'alumina_psd.nc'    ! PSD file (netCDF)
  logical, public                :: carma_emis_psd_allow_clip = .false.             ! only warn if > 1% of PSD mass is clipped
  real(r8), public               :: carma_emis_rmode          = 2.0_r8              ! lognormal mode (number median) radius (nm)
  real(r8), public               :: carma_emis_sigma          = 1.5_r8              ! lognormal geometric standard deviation
  logical, public                :: carma_alumina_rad_feedback = .true.            ! alumina in rad_climate (always in rad_diag_1,
                                                                                    ! never in rad_diag_2); acted on by build-namelist,
                                                                                    ! only reported here

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
      carma_emis_file, &
      carma_emis_scale, &
      carma_emis_psd_type, &
      carma_emis_psd_file, &
      carma_emis_psd_allow_clip, &
      carma_emis_rmode, &
      carma_emis_sigma, &
      carma_alumina_rad_feedback

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
    call mpibcast(carma_emis_file,           len(carma_emis_file),     mpichar, 0, mpicom)
    call mpibcast(carma_emis_scale,          1,                        mpir8,   0, mpicom)
    call mpibcast(carma_emis_psd_type,       len(carma_emis_psd_type), mpichar, 0, mpicom)
    call mpibcast(carma_emis_psd_file,       len(carma_emis_psd_file), mpichar, 0, mpicom)
    call mpibcast(carma_emis_psd_allow_clip, 1,                        mpilog,  0, mpicom)
    call mpibcast(carma_emis_rmode,          1,                        mpir8,   0, mpicom)
    call mpibcast(carma_emis_sigma,          1,                        mpir8,   0, mpicom)
    call mpibcast(carma_alumina_rad_feedback, 1,                       mpilog,  0, mpicom)
#endif

  end subroutine carma_model_readnl

end module carma_model_flags_mod
