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
!! Bin shift (as_do_binshift = .true.): the arsenic taken up adds to the particle
!! mass. CARMA gives every particle in bin i the mass rmass[i], so after the uptake
!! a particle of bin i has the mass
!!   m     = rmass[i] + dA_i / N_i        (dA_i = arsenic added to bin i, N_i = pc of
!!                                        the number element; the scaling cancels)
!! With j the largest bin that has rmass[j] <= m, a fraction
!!   fmove = (m - rmass[j]) / (rmass[j+1] - rmass[j])
!! of the particles goes to bin j+1 and the rest to bin j, and the arsenic core of
!! bin i (old core plus dA_i) is split with the same fractions. This conserves number,
!! total particle mass (N * rmass) and arsenic mass to round-off. DUSTnn is then smoke
!! plus arsenic, and the smoke mass is DUSTnn - DUSTASnn. Particles that would grow
!! past the largest bin stay in bin NBIN, which loses the mass above rmass[NBIN]; bin
!! NBIN has a radius of about 0.1 um and holds almost no smoke. With
!! as_do_binshift = .false. arsenic is a passive core label (the v1 behaviour).
!!
!! The module variables are set once per MPI task by CARMA_InitializeModel in
!! carma_model_mod (arsenic_msp). The routine is called once per level per CARMA step
!! from the newstate_calc.F90 override in this model directory.
!!
!! Diagnostic: the total loss rate k used at each level is kept in as_kused(iz) [1/s]
!! (0 where no rate was computed: as_gamma <= 0 or gc <= 0; the computed k is kept
!! even where k*dtime < KDT_MIN and nothing is moved). The CAM interface clears it
!! with arsenic_uptake_reset before each column and reads it with
!! arsenic_uptake_getk after CARMASTATE_Step; it is history field ASGASK. The
!! array is THREADPRIVATE because CAM may run chunks on several OpenMP threads.
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
  logical,      public :: as_do_binshift = .true.  !! add the arsenic to the particle mass (bin shift)
  real(kind=f), public :: as_gamma     = 1._f      !! uptake (sticking) coefficient, 0-1
  integer,      public :: as_igas      = 0         !! CARMA gas index of the arsenic gas
  integer,      public :: as_igroup    = 0         !! CARMA group index of the smoke
  integer,      public :: as_iecore    = 0         !! CARMA element index of the arsenic core

  ! Diagnostic: total loss rate k used at each level in the current column [1/s].
  real(kind=f), allocatable, public :: as_kused(:)
  !$OMP THREADPRIVATE(as_kused)

  public :: arsenic_uptake
  public :: arsenic_uptake_reset
  public :: arsenic_uptake_getk

contains

  !! Clear the k diagnostic (as_kused) for a new column of nz levels, allocating it
  !! if needed. Called by the CAM interface before CARMASTATE_Step for each column,
  !! so that no value from a previous column can be read back.
  subroutine arsenic_uptake_reset(nz)

    integer, intent(in) :: nz   !! number of levels in the column

    if (allocated(as_kused)) then
      if (size(as_kused) /= nz) deallocate(as_kused)
    end if
    if (.not. allocated(as_kused)) allocate(as_kused(nz))

    as_kused(:) = 0._f

    return
  end subroutine arsenic_uptake_reset


  !! Return the k diagnostic (as_kused) for the current column. Returns zeros if
  !! arsenic_uptake_reset has not been called on this thread or the sizes differ.
  subroutine arsenic_uptake_getk(nz, kused)

    integer, intent(in)       :: nz          !! number of levels in the column
    real(kind=f), intent(out) :: kused(nz)   !! total loss rate used at each level [1/s]

    kused(:) = 0._f
    if (allocated(as_kused)) then
      if (size(as_kused) == nz) kused(:) = as_kused(:)
    end if

    return
  end subroutine arsenic_uptake_getk


  !! Move arsenic from the gas (gc) to the arsenic core-mass element (pc) at level iz
  !! with an exact first-order exponential update over dtime. With as_do_binshift,
  !! then move the particles to the bins that match their new mass.
  !!
  !! NOTE: Do not put an apostrophe, or a CARMA macro name followed by a parenthesis,
  !! in comments in this file, because it is passed through the C preprocessor.
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
    real(kind=f)            :: dabin(NBIN)     ! arsenic added to each bin [scaled g/cm3]
    real(kind=f)            :: pcnum(NBIN)     ! number element after the bin shift
    real(kind=f)            :: pccore(NBIN)    ! arsenic core element after the bin shift
    real(kind=f)            :: mnew            ! particle mass after the uptake [g]
    real(kind=f)            :: fmove           ! fraction of the particles moved to bin jbin+1
    integer                 :: ibin            ! bin index
    integer                 :: jbin            ! bin that the particles of ibin move to
    integer                 :: iencn           ! number concentration element of the group

    ! Clear the k diagnostic for this level first, so that every return path below
    ! leaves 0 unless a rate is computed. Allocate it if the interface did not.
    if (.not. allocated(as_kused)) then
      allocate(as_kused(NZ))
      as_kused(:) = 0._f
    else if (size(as_kused) /= NZ) then
      deallocate(as_kused)
      allocate(as_kused(NZ))
      as_kused(:) = 0._f
    end if
    as_kused(iz) = 0._f

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

    ! Record the rate for the diagnostic, also when it is too small to act on.
    as_kused(iz) = ksum

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
      dabin(ibin) = 0._f
      if (kbinrate(ibin) > 0._f) then
        dabin(ibin) = dgc * (kbinrate(ibin) / ksum)
        pc(iz, ibin, as_iecore) = pc(iz, ibin, as_iecore) + dabin(ibin)
      end if
    end do

    if (.not. as_do_binshift) return

    ! Move the particles of each bin, with their arsenic, to the bins that match their
    ! new mass. All bins are moved from the state before the shift.
    pcnum(:)  = 0._f
    pccore(:) = 0._f

    do ibin = 1, NBIN
      if (dabin(ibin) <= 0._f) then
        pcnum(ibin)  = pcnum(ibin)  + pc(iz, ibin, iencn)
        pccore(ibin) = pccore(ibin) + pc(iz, ibin, as_iecore)
        cycle
      end if

      ! kbinrate > 0 here, so the number is > 0.
      mnew = rmass(ibin, as_igroup) + dabin(ibin) / pc(iz, ibin, iencn)

      jbin = ibin
      do while (jbin < NBIN)
        if (rmass(jbin+1, as_igroup) > mnew) exit
        jbin = jbin + 1
      end do

      if (jbin == NBIN) then
        fmove = 0._f
      else if (jbin == ibin) then
        ! The usual case. Use the added mass directly, because mnew - rmass[ibin] would
        ! cancel when the added mass is very small.
        fmove = (dabin(ibin) / pc(iz, ibin, iencn)) / &
                (rmass(ibin+1, as_igroup) - rmass(ibin, as_igroup))
      else
        fmove = (mnew - rmass(jbin, as_igroup)) / &
                (rmass(jbin+1, as_igroup) - rmass(jbin, as_igroup))
      end if

      pcnum(jbin)  = pcnum(jbin)  + (1._f - fmove) * pc(iz, ibin, iencn)
      pccore(jbin) = pccore(jbin) + (1._f - fmove) * pc(iz, ibin, as_iecore)

      if (fmove > 0._f) then
        pcnum(jbin+1)  = pcnum(jbin+1)  + fmove * pc(iz, ibin, iencn)
        pccore(jbin+1) = pccore(jbin+1) + fmove * pc(iz, ibin, as_iecore)
      end if
    end do

    pc(iz, :, iencn)     = pcnum(:)
    pc(iz, :, as_iecore) = pccore(:)

    return
  end subroutine arsenic_uptake

end module arsenic_uptake_mod
