! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! Uptake of the lumped arsenic gas (ASGAS) onto meteoric smoke particles for the
!! CARMA arsenic_msp model.
!!
!! Gas-phase arsenic collides with smoke particles at the gas-kinetic rate, with a
!! Fuchs-Sutugin transition-regime correction, and sticks with probability as_gamma.
!! The arsenic that sticks is moved into the core-mass element as_iecore of group
!! as_igroup, bin by bin in proportion to the share of each bin in the total loss rate.
!!
!! For each level iz (cgs units):
!!   n_air = rhoa / (xmet*ymet*zmet) * AVG / WTMOL_AIR                 [molec/cm3]
!!   D     = 1.76575e17 * sqrt(T) / n_air                              [cm2/s]
!!           (the kinetic-theory estimate CARMA uses for H2SO4, setupgrow.F90)
!!   vbar  = sqrt(8 RGAS T / (PI wtmol))                               [cm/s]
!!   lambda= 3 D / vbar                                                [cm]
!!   N_i   = pc(iz,ibin,ienconc) / (xmet*ymet*zmet)                    [#/cm3]
!!   r_i   = r_wet(iz,ibin,igroup)                                     [cm]
!!   Kn_i  = lambda / r_i
!!   f_i   = (1 + Kn) / (1 + (4/(3 gamma) + 0.377) Kn + 4/(3 gamma) Kn**2)
!!   k_i   = 4 PI r_i D N_i f_i                (as_do_fuchs = .true.)
!!   k_i   = gamma PI r_i**2 vbar N_i          (as_do_fuchs = .false.)
!!   k     = sum_i k_i                                                  [1/s]
!!
!!   dgc   = gc * (1 - exp(-k dtime))
!!   gc    = gc - dgc
!!   pc(iz,ibin,iecore) = pc(iz,ibin,iecore) + dgc * k_i / k
!!
!! Both gc and a core-mass pc are (mass mixing ratio * rhoa_wet), scaled by the same
!! metric factors (CARMASTATE_SetGas/GetGas and CARMASTATE_SetBin/GetBin), so moving
!! mass between them needs no unit conversion and conserves mass to round-off.
!!
!! The module variables are set once per MPI task by CARMA_InitializeModel in
!! carma_model_mod (arsenic_msp). The routine is called once per level per CARMA step
!! from the newstate_calc.F90 override in this model directory.
!!
!! @author  Asha Jain (with Claude)
!! @version Oct-2026
module arsenic_uptake_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  private

  ! Configuration, set by CARMA_InitializeModel (carma_model_mod) on every task.
  logical,      public :: as_do_uptake = .false.   !! turn the uptake on or off
  logical,      public :: as_do_fuchs  = .true.    !! apply the Fuchs-Sutugin correction
  real(kind=f), public :: as_gamma     = 1._f      !! uptake (sticking) coefficient, 0-1
  integer,      public :: as_igas      = 0         !! CARMA gas index of the arsenic gas
  integer,      public :: as_igroup    = 0         !! CARMA group index of the smoke
  integer,      public :: as_iecore    = 0         !! CARMA element index of the arsenic core

  public :: arsenic_uptake

contains

  !! Move arsenic from the gas (gc) to the arsenic core-mass element (pc) at level iz
  !! with an exact first-order exponential update over dtime.
  !!
  !! NOTE: Do not put an apostrophe in comments in this file, because it is passed
  !! through the C preprocessor.
  subroutine arsenic_uptake(carma, cstate, iz, rc)

    type(carma_type), intent(in)         :: carma   !! the carma object
    type(carmastate_type), intent(inout) :: cstate  !! the carma state object
    integer, intent(in)                  :: iz      !! z index
    integer, intent(inout)               :: rc      !! return code, negative indicates failure

    ! Below this value of k*dtime nothing is done.
    real(kind=f), parameter :: KDT_MIN  = 1.e-30_f
    ! Below this value of k*dtime, 1 - exp(-x) is evaluated with a series to avoid
    ! cancellation (relative truncation error x**3/24, i.e. < 5e-14 at the switch).
    real(kind=f), parameter :: X_SERIES = 1.e-4_f
    ! Fuchs-Sutugin constant.
    real(kind=f), parameter :: FS_C     = 0.377_f

    real(kind=f)            :: kbinrate(NBIN)  ! loss rate of gas to each bin [1/s]
    real(kind=f)            :: xyzmet          ! xmet*ymet*zmet at iz
    real(kind=f)            :: rhoa_cgs        ! air density [g/cm3]
    real(kind=f)            :: nair            ! air number density [molec/cm3]
    real(kind=f)            :: temp            ! temperature [K]
    real(kind=f)            :: diffas          ! diffusivity of the gas in air [cm2/s]
    real(kind=f)            :: vbar            ! mean thermal speed of the gas [cm/s]
    real(kind=f)            :: lamas           ! mean free path of the gas [cm]
    real(kind=f)            :: gam43           ! 4 / (3 gamma)
    real(kind=f)            :: rad             ! particle (wet) radius [cm]
    real(kind=f)            :: ndens           ! particle number density [#/cm3]
    real(kind=f)            :: knud            ! Knudsen number of the gas
    real(kind=f)            :: fcorr           ! Fuchs-Sutugin correction factor
    real(kind=f)            :: ksum            ! total loss rate [1/s]
    real(kind=f)            :: xkdt            ! ksum * dtime
    real(kind=f)            :: frac            ! 1 - exp(-xkdt)
    real(kind=f)            :: dgc             ! gas moved to particles [scaled g/cm3]
    integer                 :: ibin            ! bin index
    integer                 :: iencn           ! number concentration element of the group

    ! Check the configuration. These are cheap integer tests, and any failure is fatal.
    if ((as_igas < 1) .or. (as_igas > NGAS) .or. &
        (as_igroup < 1) .or. (as_igroup > NGROUP) .or. &
        (as_iecore < 1) .or. (as_iecore > NELEM)) then
      if (do_print) write(LUNOPRT,*) "arsenic_uptake::ERROR - bad indices: as_igas=", as_igas, &
        ", as_igroup=", as_igroup, ", as_iecore=", as_iecore
      rc = RC_ERROR
      return
    end if

    if ((itype(as_iecore) /= I_COREMASS) .or. (igelem(as_iecore) /= as_igroup)) then
      if (do_print) write(LUNOPRT,*) "arsenic_uptake::ERROR - element ", as_iecore, &
        " is not a core-mass element of group ", as_igroup
      rc = RC_ERROR
      return
    end if

    ! newstate_calc is called twice per step (cloudy and clear) when do_incloud is
    ! set, which would apply the uptake twice. carma_model_mod refuses that setting
    ! at initialization; this is a second guard.
    if (do_incloud) then
      if (do_print) write(LUNOPRT,*) "arsenic_uptake::ERROR - do_incloud is not supported."
      rc = RC_ERROR
      return
    end if

    ! Early returns.
    if (as_gamma <= 0._f) return
    if (gc(iz, as_igas) <= 0._f) return

    ! Gas properties.
    xyzmet   = xmet(iz) * ymet(iz) * zmet(iz)
    rhoa_cgs = rhoa(iz) / xyzmet
    nair     = rhoa_cgs * AVG / WTMOL_AIR
    temp     = t(iz)
    diffas   = 1.76575e+17_f * sqrt(temp) / nair
    vbar     = sqrt(8._f * RGAS * temp / (PI * gwtmol(as_igas)))
    lamas    = 3._f * diffas / vbar
    gam43    = 4._f / (3._f * as_gamma)

    ! Loss rate to each bin.
    iencn = ienconc(as_igroup)
    ksum  = 0._f

    do ibin = 1, NBIN
      kbinrate(ibin) = 0._f

      rad   = r_wet(iz, ibin, as_igroup)
      ndens = pc(iz, ibin, iencn) / xyzmet

      if ((rad <= 0._f) .or. (ndens <= 0._f)) cycle

      if (as_do_fuchs) then
        knud  = lamas / rad
        fcorr = (1._f + knud) / (1._f + (gam43 + FS_C) * knud + gam43 * knud * knud)
        kbinrate(ibin) = 4._f * PI * rad * diffas * ndens * fcorr
      else
        kbinrate(ibin) = as_gamma * PI * rad * rad * vbar * ndens
      end if

      ksum = ksum + kbinrate(ibin)
    end do

    xkdt = ksum * dtime
    if (xkdt < KDT_MIN) return

    ! Fraction of the gas taken up during the step, 1 - exp(-xkdt).
    if (xkdt < X_SERIES) then
      frac = xkdt * (1._f - xkdt * (0.5_f - xkdt / 6._f))
    else
      frac = 1._f - exp(-xkdt)
    end if

    ! Move the mass from the gas to the arsenic core of each bin.
    dgc = gc(iz, as_igas) * frac
    gc(iz, as_igas) = gc(iz, as_igas) - dgc

    do ibin = 1, NBIN
      if (kbinrate(ibin) > 0._f) then
        pc(iz, ibin, as_iecore) = pc(iz, ibin, as_iecore) + dgc * (kbinrate(ibin) / ksum)
      end if
    end do

    return
  end subroutine arsenic_uptake

end module arsenic_uptake_mod
