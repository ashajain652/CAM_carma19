!! Minimal stand-ins for the CAM modules that CAM's CARMA wrapper modules
!! (cam/carma_precision_mod.F90, cam/carma_constants_mod.F90) use, so CARMA's
!! own code can be compiled outside of CAM. Values are the CESM ones:
!!   physconst            -> shr_const_mod (as in CAM's physconst.F90)
!!   radconstants         -> band limits from src/physics/rrtmg/radconstants.F90
!!   cam_history_support  -> fillvalue from src/control/cam_history_support.F90

module physconst
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod
  implicit none
  real(r8), public, parameter :: avogad      = shr_const_avogad
  real(r8), public, parameter :: boltz       = shr_const_boltz
  real(r8), public, parameter :: latice      = shr_const_latice
  real(r8), public, parameter :: latvap      = shr_const_latvap
  real(r8), public, parameter :: pi          = shr_const_pi
  real(r8), public, parameter :: r_universal = shr_const_rgas
  real(r8), public, parameter :: rhoh2o      = shr_const_rhofw
end module physconst

module radconstants
  use shr_kind_mod,  only: r8 => shr_kind_r8
  implicit none
  integer, public, parameter :: nswbands = 14
  integer, public, parameter :: nlwbands = 16

  real(r8), parameter :: wavenum_low(nswbands) = & ! in cm^-1
    (/2600._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, &
      8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,  820._r8/)
  real(r8), parameter :: wavenum_high(nswbands) = & ! in cm^-1
    (/3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, 8050._r8, &
     12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8, 2600._r8/)
  real(r8), parameter :: wavenumber1_longwave(nlwbands) = &
    (/   10._r8,  350._r8, 500._r8,   630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, &
       1180._r8, 1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2600._r8 /)
  real(r8), parameter :: wavenumber2_longwave(nlwbands) = &
    (/  350._r8,  500._r8,  630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, 1180._r8, &
       1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2600._r8, 3250._r8 /)

contains

  ! Same as radconstants.F90 with units = 'cm'.
  subroutine get_lw_spectral_boundaries_cm(low_boundaries, high_boundaries)
    real(r8), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
    low_boundaries  = 1._r8/wavenumber2_longwave
    high_boundaries = 1._r8/wavenumber1_longwave
  end subroutine get_lw_spectral_boundaries_cm

  subroutine get_sw_spectral_boundaries_cm(low_boundaries, high_boundaries)
    real(r8), intent(out) :: low_boundaries(nswbands), high_boundaries(nswbands)
    low_boundaries  = 1._r8/wavenum_high
    high_boundaries = 1._r8/wavenum_low
  end subroutine get_sw_spectral_boundaries_cm
end module radconstants

module cam_history_support
  use shr_kind_mod,  only: r8 => shr_kind_r8
  implicit none
  real(r8), parameter, public :: fillvalue = 1.e36_r8
end module cam_history_support
