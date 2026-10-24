!! This CARMA model is for reentry arsenic carried on meteor smoke aerosols. The
!! meteor smoke part is the meteor_smoke model, which is based upon Bardeen et al.,
!! JGR, 2008.
!!
!! In addition to the smoke, there is:
!!
!!   - a lumped arsenic gas (ASGAS), counted as mass of As atoms, that is emitted
!!     uniformly in time and per unit area into the model layer containing
!!     carma_as_emis_alt (km above sea level)
!!   - a core mass element (DUSTAS) in the smoke group, that holds the arsenic that
!!     has been taken up by the smoke. The uptake itself is done on the CARMA side
!!     by arsenic_uptake_mod, which is configured in CARMA_InitializeModel().
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitGas()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeModel()
!!
!! NOTE: CARMA_EmitGas() is not called by the standard carma_intr.F90. It is called
!! from the override copy of carma_intr.F90 in this model directory.
!!
!! @version Jan-2011
!! @author  Chuck Bardeen
!! @version Oct-2026 (arsenic_msp, based upon meteor_smoke)
module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use carma_flags_mod
  use carma_model_flags_mod

  use spmd_utils,     only: masterproc
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radconstants,   only: nswbands, nlwbands
  use cam_abortutils, only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc

#if ( defined SPMD )
  use mpishorthand
#endif  

  implicit none

  private

  ! Declare the public methods.
  public CARMA_DefineModel
  public CARMA_Detrain
  public CARMA_DiagnoseBins
  public CARMA_DiagnoseBulk
  public CARMA_EmitGas
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_WetDeposition

  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 1               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 2               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 28              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 1               !! Number of gases

  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)

  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_METEOR_SMOKE   = 1       !! meteor smoke
  integer, public, parameter      :: I_ARSENIC        = 2       !! arsenic (on smoke)

  ! Define any gas compositions that are used.
  !
  ! NOTE: This is outside of the range of the gas compositions that CARMA treats
  ! specially (I_GCOMP_H2O, I_GCOMP_H2SO4, I_GCOMP_SO2), so CARMA does nothing special
  ! with it.
  integer, public, parameter      :: I_GCOMP_AS       = 101     !! arsenic (lumped gas)

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_DUST     = 1         !! meteor smoke

  integer, public, parameter      :: I_ELEM_DUST    = 1         !! meteor smoke
  integer, public, parameter      :: I_ELEM_ASCORE  = 2         !! arsenic on meteor smoke (core mass)

  integer, public, parameter      :: I_GAS_AS       = 1         !! arsenic (lumped gas)

  ! Densities and molecular weights.
  real(kind=f), public, parameter :: RHO_METEOR_SMOKE = 2.0_f      !! density of meteor smoke particles (g/cm3)
  real(kind=f), public, parameter :: RHO_ARSENIC      = 5.73_f     !! density of arsenic (g/cm3)
  real(kind=f), public, parameter :: WTMOL_AS         = 74.9216_f  !! molecular weight of arsenic gas, as As atoms (g/mol)

  ! These variables are all set during initialization and are used to calculate
  ! emission tendencies.
  integer                             :: carma_emis_nLevs        ! number of emission levels
  real(r8), allocatable, dimension(:) :: carma_emis_lev          ! emission levels (Pa)
  real(r8), allocatable, dimension(:) :: carma_emis_rate         ! emission rate lookup table (# cm-3 s-1)
  integer                             :: carma_emis_ilev_min     ! index of minimum level in table
  integer                             :: carma_emis_ilev_max     ! index of maximum level in table
  integer                             :: carma_emis_ilev_incr    ! index increment to increase level
  real(r8)                            :: carma_emis_expected     ! Expected emission rate per column (kg/m2/s)

  real(r8)                            :: carma_as_flux = 0._r8   ! Arsenic gas emission flux per column (kg/m2/s)

contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DefineModel(carma, rc)
    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    real(kind=f), parameter            :: rmin     = 2e-8_f   ! minimum radius (cm)
    real(kind=f), parameter            :: vmrat    = 2.0_f    ! volume ratio

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?
    
    ! Default return code.
    rc = RC_OK
    
    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_emis_total    = ', carma_emis_total
      if (do_print) write(LUNOPRT,*) '  carma_emis_file     = ', trim(carma_emis_file)
      if (do_print) write(LUNOPRT,*) '  carma_as_emis_total = ', carma_as_emis_total, ' (Gg/yr)'
      if (do_print) write(LUNOPRT,*) '  carma_as_emis_alt   = ', carma_as_emis_alt, ' (km)'
      if (do_print) write(LUNOPRT,*) '  carma_as_gamma      = ', carma_as_gamma
      if (do_print) write(LUNOPRT,*) '  carma_as_do_uptake  = ', carma_as_do_uptake
      if (do_print) write(LUNOPRT,*) '  carma_as_do_fuchs   = ', carma_as_do_fuchs
      if (do_print) write(LUNOPRT,*) '  carma_as_do_binshift= ', carma_as_do_binshift
      if (do_print) write(LUNOPRT,*) '  carma_as_top_noflux = ', carma_as_top_noflux
    end if
    
    
    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    call CARMAGROUP_Create(carma, I_GRP_DUST, "meteor smoke", rmin, vmrat, I_SPHERE, 1._f, .false., &
                          rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, shortname="DUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_DUST, I_GRP_DUST, "meteor smoke", RHO_METEOR_SMOKE, &
         I_INVOLATILE, I_METEOR_SMOKE, rc, shortname="DUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    ! The arsenic taken up by the smoke is a core mass in the smoke group. It moves
    ! with the smoke (coagulation, sedimentation, deposition), but it does not change
    ! the size of the particles.
    call CARMAELEMENT_Create(carma, I_ELEM_ASCORE, I_GRP_DUST, "arsenic on smoke", RHO_ARSENIC, &
         I_COREMASS, I_ARSENIC, rc, shortname="DUSTAS")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')


    ! Define the Solutes


    ! Define the Gases
    !
    ! NOTE: CARMA does not accept an ivaprtn of 0. Even with do_grow false, supersat()
    ! (and so vaporp()) is called for every gas in CARMASTATE_CreateFromReference and
    ! CARMASTATE_Step, and vaporp() returns RC_ERROR for an unknown vapor pressure
    ! routine. Murphy & Koop [2005] is used as a placeholder, since it only depends
    ! upon temperature, never fails and does not need a water vapor gas (unlike
    ! H2SO4 Ayers [1980]). The resulting arsenic supersaturations are meaningless, but
    ! they are only used by growth, which is not allowed for this model (see
    ! CARMA_InitializeModel).
    call CARMAGAS_Create(carma, I_GAS_AS, "arsenic (lumped gas)", WTMOL_AS, I_VAPRTN_H2O_MURPHY2005, &
         I_GCOMP_AS, rc, shortname="ASGAS")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')


    ! Define the Processes
    call CARMA_AddCoagulation(carma, I_GRP_DUST, I_GRP_DUST, I_GRP_DUST, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')
    
    return
  end subroutine CARMA_DefineModel


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!
  !!  @see CARMASTATE_SetDetrain
  subroutine CARMA_Detrain(carma, cstate, cam_in, dlf, state, icol, dt, rc, rliq, prec_str, snow_str, &
     tnd_qsnow, tnd_nsnow)
    use camsrfexch,         only: cam_in_t
    use physconst,          only: latice, latvap, cpair

    implicit none

    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in           !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(r8), intent(in)                 :: dt               !! time step (s)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional      :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(out), optional      :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)

    ! Default return code.
    rc = RC_OK
    
    return
  end subroutine CARMA_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, rliq, prec_str, snow_str)
    use time_manager,     only: is_first_step

    implicit none

    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)       :: state        !! physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    
    real(r8)                             :: mmr(pver) !! elements mass mixing ratio
    integer                              :: ibin      !! bin index
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the mass in each bin from the CAM state.
    
    return
  end subroutine CARMA_DiagnoseBins


  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    type(cam_out_t),      intent(inout)  :: cam_out   !! cam output to surface models
    type(physics_state), intent(in)      :: state     !! physics state variables
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer
    type(physics_ptend), intent(inout)   :: ptend     !! constituent tendencies
    integer, intent(in)                  :: icol      !! column index
    real(r8), intent(in)                 :: dt        !! time step
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(inout), optional    :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(inout), optional    :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    real(r8), intent(out), optional      :: re_ice(pcols,pver)    !! ice effective radius (m)
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the bulk mass from the CARMA state.
    
    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA gases. This is called from the override copy
  !! of carma_intr.F90 in this model directory.
  !!
  !! The arsenic gas source is constant in time and the same per unit area everywhere.
  !! In each column, the whole flux carma_as_flux (kg/m2/s) goes into the single layer
  !! that contains carma_as_emis_alt (km above sea level). The height above sea level
  !! is zi + phis/gravit, i.e. geopotential height.
  !!
  !! @author  arsenic_msp, based upon CARMA_EmitParticle by Chuck Bardeen
  !! @version Oct-2026
  subroutine CARMA_EmitGas(carma, igas, icnst, dt, state, cam_in, tendency, rc)
    use camsrfexch,    only: cam_in_t
    use physconst,     only: gravit
    use cam_logfile,   only: iulog

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: igas                  !! gas index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ncol                  ! number of columns in chunk
    integer                            :: icol                  ! column index
    integer                            :: k                     ! vertical index
    integer                            :: nmiss                 ! number of columns with no emission layer
    integer                            :: icol_miss             ! first column with no emission layer
    logical                            :: found                 ! was the emission layer found?
    real(r8)                           :: zasl                  ! emission height above the surface (m)
    real(r8)                           :: zasl_miss             ! emission height above the surface in icol_miss (m)

    ! Default return code.
    rc = RC_OK

    ncol = state%ncol

    ! For emissions into the atmosphere, put the emission here.
    tendency(:, :) = 0.0_r8

    if (igas == I_GAS_AS) then

      nmiss     = 0
      icol_miss = 0
      zasl_miss = 0.0_r8

      do icol = 1, ncol

        ! The emission altitude is above sea level, but zi is above the surface.
        zasl  = carma_as_emis_alt * 1000._r8 - state%phis(icol) / gravit
        found = .false.

        ! Put the whole column flux into the one layer containing zasl. The interval
        ! is closed at the bottom and open at the top, so only one layer can match.
        !
        ! NOTE: The tendency is on a mass mixing ratio (kg/kg/s).
        do k = 1, pver
          if ((state%zi(icol, k+1) <= zasl) .and. (zasl < state%zi(icol, k))) then
            tendency(icol, k) = carma_as_flux / (state%pdel(icol, k) / gravit)
            found = .true.
            exit
          end if
        end do

        if (.not. found) then
          nmiss = nmiss + 1

          if (nmiss == 1) then
            icol_miss = icol
            zasl_miss = zasl
          end if
        end if
      end do

      ! Every column should have an emission layer, since CARMA_InitializeModel checks
      ! that the emission altitude is below the model top.
      if (nmiss > 0) then
        write(iulog,*) 'CARMA_EmitGas: ERROR - no emission layer found in ', nmiss, ' of ', ncol, &
          ' columns in chunk ', state%lchnk
        write(iulog,*) 'CARMA_EmitGas:   first column = ', icol_miss, ', lat = ', state%lat(icol_miss) * RAD2DEG, &
          ', lon = ', state%lon(icol_miss) * RAD2DEG
        write(iulog,*) 'CARMA_EmitGas:   emission height above surface (m) = ', zasl_miss, &
          ', model top (m) = ', state%zi(icol_miss, 1), ', surface (m) = ', state%zi(icol_miss, pver+1)
        call endrun('CARMA_EmitGas: no layer contains carma_as_emis_alt.')
      end if
    end if

    return
  end subroutine CARMA_EmitGas


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Chuck Bardeen
  !! @version Jan-2011
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use camsrfexch,       only: cam_in_t
    use physconst,     only: gravit

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ncol                  ! number of columns in chunk
    integer                            :: icol                  ! column index
    integer                            :: igroup                ! the index of the carma aerosol group
    integer                            :: k                     ! vertical index
    integer                            :: ilev                  ! level index in emissions data
    real(r8)                           :: r(NBIN)               ! bin center
    real(r8)                           :: dr(NBIN)              ! bin width
    real(r8)                           :: rmass(NBIN)           ! bin mass
    real(r8)                           :: pressure              ! pressure (Pa)
    real(r8)                           :: thickness             ! layer thickness (m)
    real(r8)                           :: rate                  ! emission rate (#/cm-3/s)
    real(r8)                           :: massflux              ! emission mass flux (kg/m2/s)
    real(r8)                           :: columnMass            ! mass of the total column (kg/m2/s)
    real(r8)                           :: scale                 ! scaling factor to conserve the expected mass


    ! Default return code.
    rc = RC_OK

    ncol = state%ncol

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.            
    tendency(:ncol, :pver) = 0.0_r8


    ! Only do emission for the first bin of the meteor smoke element.
    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup)
    if (RC < RC_ERROR) return

    call CARMAGROUP_GET(carma, igroup, rc, r=r, dr=dr, rmass=rmass)
    if (RC < RC_ERROR) return

    ! For meteoritic dust, the source from the smoke only goes into the
    ! smallest bin (~1.3 nm). The depth that the micrometeorite penetrates
    ! is proportional to the pressure, so the emission is a function of
    ! pressure.
    !
    ! NOTE: This tests the element rather than the group shortname, since the
    ! arsenic core element (DUSTAS) is also in the DUST group. Testing the group
    ! would put the smoke mass into DUSTAS01 too.
    if ((ielem == I_ELEM_DUST) .and. (ibin == 1)) then

      ! Set tendencies for any sources or sinks in the atmosphere.
      do k = 1, pver
        do icol = 1, ncol
      
          pressure = state%pmid(icol, k)
          
          ! This is roughly a log-normal approximation to the production
          ! rate, but only applies from about 70 to 110 km.
          !
          ! NOTE: Based upon US Standard Atmosphere 1976.
          if ((pressure >= carma_emis_lev(carma_emis_ilev_min)) .and. &
              (pressure <= carma_emis_lev(carma_emis_ilev_max))) then

            ! The rates are in terms of # cm-3 s-1, but were really derived
            ! from the mass flux of meteoritic dust. Since we are using a
            ! size different that 1.3 nm for the smallest bin, scale the
            ! number appropriately.
            !
            ! The values are in a lookup table, so find the two numbers
            ! surrounding the pressure and do a linear interpolation on the
            ! rate. This linear search is kind of expensive, particularly if
            ! there are a lot of points.
            ! 
            ! NOTE: The tendency is on a mass mixing ratio (kg/kg/s)
            do ilev = carma_emis_ilev_min, (carma_emis_ilev_max - carma_emis_ilev_incr), carma_emis_ilev_incr
              if ((pressure >= carma_emis_lev(ilev)) .and. (pressure <= carma_emis_lev(ilev+carma_emis_ilev_incr))) then
                rate = carma_emis_rate(ilev)
                
                if (pressure > carma_emis_lev(ilev)) then
                  rate = rate + &
                    ((carma_emis_rate(ilev+carma_emis_ilev_incr) - &
                    carma_emis_rate(ilev)) / (carma_emis_lev(ilev+carma_emis_ilev_incr) - &
                    carma_emis_lev(ilev))) * (pressure - carma_emis_lev(ilev))
                end if
                
                rate = rate * (((1.3e-7_r8)**3) / (r(ibin)**3))
                exit
              end if
            end do
            
            ! Calculate the mass flux in terms of kg/m3/s
            massflux = (rate * rmass(ibin) * 1.0e-3_r8 * 1.0e6_r8)
            
            ! Convert the mass flux to a tendency on the mass mixing ratio.
            thickness = state%zi(icol, k) - state%zi(icol, k+1)
            tendency(icol, k) = (massflux * thickness) / (state%pdel(icol, k) / gravit)        
          end if
        enddo
      enddo

      ! Scale the columns to keep the total mass influx in the column a
      ! constant.
      !
      ! NOTE: Guard against a column with no emission levels, which would otherwise
      ! give a divide by zero.
      do icol = 1, ncol
        columnMass = sum(tendency(icol, :) * (state%pdel(icol, :) / gravit))
        
        if (columnMass > 0._r8) then
          scale = carma_emis_expected / columnMass
          tendency(icol, :) = tendency(icol, :) * scale
        else
          tendency(icol, :) = 0.0_r8
        end if
      end do
    end if
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use ioFileMod,    only: getfil
    use constituents, only: pcnst
    use wrap_nf
    use physconst,    only: p_rearth=>rearth
    use ppgrid,       only: pverp
    use ref_pres,     only: pref_edge
    use arsenic_uptake_mod, only: as_do_uptake, as_do_fuchs, as_do_binshift, as_gamma, as_igas, as_igroup, as_iecore

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer                            :: ilev                  ! level index
    integer                            :: fid                   ! file id
    integer                            :: lev_did               ! level dimension id
    integer                            :: lev_vid               ! level variable id
    integer                            :: rate_vid              ! rate variable
    integer                            :: tmp
    character(len=256)                 :: efile                 ! emission file name
    real(r8)                           :: ztop                  ! approximate height of the model top (m)

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?

    ! Default return code.
    rc = RC_OK

    ! Add initialization here.
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")
    
    ! Initialize the emissions rate table.
    if (carma_do_emission) then 
      if (masterproc) then

        ! Open the netcdf file (read only)
        call getfil(carma_emis_file, efile, fid)
        if (do_print) write(LUNOPRT,*) 'carma_init(): Reading particle emission rates from ', efile

        call wrap_open(efile, 0, fid)

        ! Alocate the table arrays
        call wrap_inq_dimid(fid, "lev", lev_did)
        call wrap_inq_dimlen(fid, lev_did, carma_emis_nLevs)
      endif
    
#if ( defined SPMD )
      call mpibcast(carma_emis_nLevs, 1, mpiint, 0, mpicom)
#endif

      allocate(carma_emis_lev(carma_emis_nLevs))
      allocate(carma_emis_rate(carma_emis_nLevs))

      if (masterproc) then
        ! Read in the tables.
        call wrap_inq_varid(fid, 'MSMOKE', rate_vid)
        call wrap_get_var_realx(fid, rate_vid, carma_emis_rate)

        call wrap_inq_varid(fid, 'lev', lev_vid)
        call wrap_get_var_realx(fid, lev_vid, carma_emis_lev)

        ! Close the file.
        call wrap_close(fid)

        ! Find out where the bounds of the table are and in what order
        ! the pressures levels are in.
        carma_emis_ilev_min = 1
        carma_emis_ilev_max = carma_emis_nLevs

        do ilev = 1, carma_emis_nLevs
          if (carma_emis_rate(ilev) <= 0.0) then
            carma_emis_ilev_min  = ilev + 1
          else
            exit  
          endif
        end do

        do ilev = carma_emis_nLevs, 1, -1
          if (carma_emis_rate(ilev) <= 0.0) then
            carma_emis_ilev_max  = ilev - 1
          else
            exit  
          endif
        end do

        if (carma_emis_lev(carma_emis_ilev_min) < carma_emis_lev(carma_emis_ilev_max)) then
          carma_emis_ilev_incr = 1
        else
          carma_emis_ilev_incr = -1
          tmp = carma_emis_ilev_min
          carma_emis_ilev_min = carma_emis_ilev_max
          carma_emis_iLev_max = tmp 
        endif

        if (do_print) write(LUNOPRT,*) ''
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_nLevs     = ', carma_emis_nLevs
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_min  = ', carma_emis_ilev_min 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_max  = ', carma_emis_ilev_max 
        if (do_print) write(LUNOPRT,*) 'carma_init(): carma_emis_ilev_incr = ', carma_emis_ilev_incr 
        if (do_print) write(LUNOPRT,*) ''
        
        if (do_print) write(LUNOPRT,*) 'level, pressure (Pa), emission rate (# cm-3 sec-1)'
        do ilev = carma_emis_ilev_min, carma_emis_ilev_max, carma_emis_ilev_incr
          if (do_print) write(LUNOPRT,*) ilev, carma_emis_lev(ilev), carma_emis_rate(ilev)
        enddo
        
        if (do_print) write(LUNOPRT, *) 'carma_init(): Total Emission = ', carma_emis_total, ' (kt/yr)'
        carma_emis_expected = ((carma_emis_total * 1e6_r8) / (3600.0_r8 * 24.0_r8 * 365.0_r8)) / &
             (4.0_r8 * PI * ((REARTH / 100._r8) ** 2))
        if (do_print) write(LUNOPRT,*) 'carma_init(): Done with emission table.'

      endif

#if ( defined SPMD )
      call mpibcast(carma_emis_lev,  carma_emis_nLevs, mpir8, 0, mpicom)
      call mpibcast(carma_emis_rate, carma_emis_nLevs, mpir8, 0, mpicom)

      call mpibcast(carma_emis_expected,  1, mpir8,  0, mpicom)

      call mpibcast(carma_emis_ilev_min,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_max,  1, mpiint, 0, mpicom)
      call mpibcast(carma_emis_ilev_incr, 1, mpiint, 0, mpicom)
#endif

    endif
    

    ! Check that the CARMA configuration is one that the arsenic uptake supports.
    !
    ! NOTE: The uptake is done in newstate_calc after coagulation, once per
    ! newstate_calc call. Growth is not supported (the arsenic gas has a placeholder
    ! vapor pressure routine and the uptake is not part of the growth substepping), and
    ! with in-cloud calculations newstate_calc is called more than once per step.
    if (carma_do_grow) then
      call endrun("CARMA_InitializeModel: arsenic_msp requires carma_do_grow=.false.")
    end if

    if (carma_do_incloud) then
      call endrun("CARMA_InitializeModel: arsenic_msp requires carma_do_incloud=.false.")
    end if

    if ((carma_as_gamma < 0._r8) .or. (carma_as_gamma > 1._r8)) then
      call endrun("CARMA_InitializeModel: carma_as_gamma must be between 0 and 1.")
    end if


    ! Initialize the arsenic gas emission.
    !
    ! The source is constant in time and the same per unit area everywhere, so the
    ! flux per column is the global total (Gg/yr -> kg/s, noleap year) divided by the
    ! area of the earth.
    carma_as_flux = ((carma_as_emis_total * 1e6_r8) / (3600.0_r8 * 24.0_r8 * 365.0_r8)) / &
         (4.0_r8 * PI * (p_rearth ** 2))

    if (masterproc) then
    
      ! As a rough guide, make sure that the emission altitude is below the model top.
      ! This uses a 7 km scale height with the reference pressures, which underestimates
      ! the height of the top of a model that extends into the thermosphere. The layer
      ! is found exactly in each column in CARMA_EmitGas.
      ztop = 7000._r8 * log(pref_edge(pverp) / pref_edge(1))

      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'carma_init(): Arsenic gas emission:'
      if (do_print) write(LUNOPRT,*) 'carma_init():   carma_as_emis_total = ', carma_as_emis_total, ' (Gg/yr)'
      if (do_print) write(LUNOPRT,*) 'carma_init():   carma_as_emis_alt   = ', carma_as_emis_alt, ' (km)'
      if (do_print) write(LUNOPRT,*) 'carma_init():   carma_as_flux       = ', carma_as_flux, ' (kg/m2/s)'
      if (do_print) write(LUNOPRT,*) 'carma_init():   approx. model top   = ', ztop / 1000._r8, ' (km)'

      if ((carma_as_emis_alt * 1000._r8) >= ztop) then
        call endrun("CARMA_InitializeModel: carma_as_emis_alt is at or above the model top.")
      end if
    end if


    ! Configure the arsenic uptake, which is done on the CARMA side. This needs to be
    ! done on every task.
    as_do_uptake = carma_as_do_uptake
    as_do_fuchs  = carma_as_do_fuchs
    as_do_binshift = carma_as_do_binshift
    as_gamma     = carma_as_gamma
    as_igas      = I_GAS_AS
    as_igroup    = I_GRP_DUST
    as_iecore    = I_ELEM_ASCORE

    if (masterproc) then
      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'carma_init(): Arsenic uptake:'
      if (do_print) write(LUNOPRT,*) 'carma_init():   as_do_uptake = ', as_do_uptake
      if (do_print) write(LUNOPRT,*) 'carma_init():   as_do_fuchs  = ', as_do_fuchs
      if (do_print) write(LUNOPRT,*) 'carma_init():   as_do_binshift = ', as_do_binshift
      if (do_print) write(LUNOPRT,*) 'carma_init():   as_gamma     = ', as_gamma
      if (do_print) write(LUNOPRT,*) 'carma_init():   as_igas      = ', as_igas
      if (do_print) write(LUNOPRT,*) 'carma_init():   as_igroup    = ', as_igroup
      if (do_print) write(LUNOPRT,*) 'carma_init():   as_iecore    = ', as_iecore
      if (do_print) write(LUNOPRT,*) ''
    end if
    
    return
  end subroutine CARMA_InitializeModel


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeParticle(carma, ielem, ibin, latvals, lonvals, mask, q, rc)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plev, plon

    implicit none

    type(carma_type), intent(in)  :: carma      !! the carma object
    integer,          intent(in)  :: ielem      !! element index
    integer,          intent(in)  :: ibin       !! bin index
    real(r8),         intent(in)  :: latvals(:) !! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) !! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    !! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     !! mass mixing ratio (gcol, lev)
    integer,          intent(out) :: rc         !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    !
    ! NOTE: Initialized to 0. by the caller, so nothing needs to be done.

    return
  end subroutine CARMA_InitializeParticle
  
  
  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure
    
    integer    :: icol
 
    ! Default return code.
    rc = RC_OK
    
    return
  end subroutine CARMA_WetDeposition 

end module
