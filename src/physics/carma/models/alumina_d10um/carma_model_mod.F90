!! This CARMA model is for involatile alumina (Al2O3) particles from spacecraft
!! reentry ablation. It is based upon the meteor_smoke model.
!!
!! PARTICLE SIZE SENSITIVITY MODEL alumina_d10um: a copy of models/alumina with rmin set to
!! the monodisperse emission radius (5 um, i.e. 10 um diameter) and 16 bins, so all
!! emitted mass goes into bin 1. Only NBIN and rmin differ from models/alumina; the
!! other modules are symlinks to ../alumina. Default PSD file (build-namelist):
!! atm/waccm/emis/alumina_mono_d10um_psd.nc.
!!
!! Emissions are discrete point events read from a netCDF event file (see
!! alumina_events_mod). All of an event row's mass is injected within one model
!! timestep into the nearest column and the layer containing the row's altitude,
!! and it is split across the CARMA bins using an altitude dependent particle size
!! distribution (see alumina_psd_mod).
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!   - CARMA_InitializeModel()
!!
!! @version Sep-2026
!! @author  Chuck Bardeen (meteor_smoke template)
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
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_WetDeposition

  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 1               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 1               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 16              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 0               !! Number of gases

  ! Relative humidities for the mie calculations. Alumina does not swell, so only the
  ! first value is used when the optics files are created.
  integer, public, parameter      :: NMIE_RH  = 8               !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH) = (/ 0._f, 0.5_f, 0.7_f, 0.8_f, 0.9_f, 0.95_f, 0.98_f, 0.99_f /)

  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_ALUMINA      = 1         !! alumina (Al2O3)

  ! Define group, element, solute and gas indexes.
  integer, public, parameter      :: I_GRP_ALUM     = 1         !! alumina

  integer, public, parameter      :: I_ELEM_ALUM    = 1         !! alumina

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
    real(kind=f), parameter            :: RHO_ALUMINA = 3.95_f  ! density of alumina particles (g/cm3)
    real(kind=f), parameter            :: rmin     = 5.0000e-04_f   ! minimum radius (cm) = emission radius 5 um
    real(kind=f), parameter            :: vmrat    = 2.0_f    ! volume ratio

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?
    complex(kind=f)                    :: refidx(NWAVE)         ! refractive indices

    ! Default return code.
    rc = RC_OK

    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_DefineModel: CARMA_Get failed.")

      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_emis_file           = ', trim(carma_emis_file)
      if (do_print) write(LUNOPRT,*) '  carma_emis_scale          = ', carma_emis_scale
      if (do_print) write(LUNOPRT,*) '  carma_emis_psd_type       = ', trim(carma_emis_psd_type)
      if (do_print) write(LUNOPRT,*) '  carma_emis_psd_file       = ', trim(carma_emis_psd_file)
      if (do_print) write(LUNOPRT,*) '  carma_emis_psd_allow_clip = ', carma_emis_psd_allow_clip
      if (do_print) write(LUNOPRT,*) '  carma_emis_rmode (nm)     = ', carma_emis_rmode
      if (do_print) write(LUNOPRT,*) '  carma_emis_sigma          = ', carma_emis_sigma
      if (do_print) write(LUNOPRT,*) '  carma_alumina_rad_feedback= ', carma_alumina_rad_feedback
    end if


    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    !
    ! Al2O3 refractive indices at the CARMA band centres (the 16 RRTMG LW bands, then
    ! the 14 SW bands), interpolated from tropf_real_alumina_2.csv (n) and
    ! tropf_kapoor_combined_img_alumina.csv (k). Generated by
    ! ~/MIT/cesm2_port/emissions/make_alumina_refidx.py; regenerate rather than edit.
    if (NWAVE /= 30) call endrun('CARMA_DefineModel: alumina refidx table needs NWAVE = 30.')

    refidx(:) = (/ &
      (3.009408_f, 1.711257e-03_f), &   ! LW01  514.2857 um
      (4.946195_f, 2.490638e-01_f), &   ! LW02   24.2857 um
      (7.520094_f, 9.899022e-01_f), &   ! LW03   17.9365 um
      (0.408993_f, 3.054029e+00_f), &   ! LW04   15.0794 um
      (0.074378_f, 1.770150e+00_f), &   ! LW05   13.2404 um
      (0.148085_f, 4.233715e-01_f), &   ! LW06   11.1996 um
      (0.921320_f, 4.577217e-02_f), &   ! LW07    9.7317 um
      (1.147377_f, 2.209648e-02_f), &   ! LW08    8.8669 um
      (1.314367_f, 8.386046e-03_f), &   ! LW09    7.8344 um
      (1.432675_f, 3.379822e-03_f), &   ! LW10    6.9755 um
      (1.502389_f, 1.510870e-03_f), &   ! LW11    6.1562 um
      (1.579544_f, 8.161219e-04_f), &   ! LW12    5.1816 um
      (1.605277_f, 7.373231e-04_f), &   ! LW13    4.6261 um
      (1.621110_f, 6.931451e-04_f), &   ! LW14    4.3143 um
      (1.637415_f, 6.533971e-04_f), &   ! LW15    4.0151 um
      (1.657288_f, 5.897629e-04_f), &   ! LW16    3.4615 um
      (1.657288_f, 5.897629e-04_f), &   ! SW01    3.4615 um
      (1.681661_f, 5.177378e-04_f), &   ! SW02    2.7885 um
      (1.691818_f, 4.740957e-04_f), &   ! SW03    2.3253 um
      (1.698969_f, 4.521094e-04_f), &   ! SW04    2.0461 um
      (1.706502_f, 4.216457e-04_f), &   ! SW05    1.7839 um
      (1.717394_f, 3.782490e-04_f), &   ! SW06    1.4624 um
      (1.719662_f, 3.782490e-04_f), &   ! SW07    1.2705 um
      (1.728741_f, 3.782490e-04_f), &   ! SW08    1.0102 um
      (1.739828_f, 4.290628e-04_f), &   ! SW09    0.7016 um
      (1.739828_f, 5.396475e-04_f), &   ! SW10    0.5333 um
      (1.752836_f, 8.542263e-04_f), &   ! SW11    0.3932 um
      (1.781104_f, 1.503466e-03_f), &   ! SW12    0.3040 um
      (1.845483_f, 3.702449e-03_f), &   ! SW13    0.2316 um
      (1.289398_f, 1.007862e-02_f) /)   ! SW14    8.0206 um

    call CARMAGROUP_Create(carma, I_GRP_ALUM, "alumina", rmin, vmrat, I_SPHERE, 1._f, .false., &
                          rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, shortname="CRALUM", refidx=refidx, do_mie=.true.)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')


    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_ALUM, I_GRP_ALUM, "alumina", RHO_ALUMINA, &
         I_INVOLATILE, I_ALUMINA, rc, shortname="CRALUM")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')


    ! Define the Solutes


    ! Define the Gases


    ! Define the Processes
    call CARMA_AddCoagulation(carma, I_GRP_ALUM, I_GRP_ALUM, I_GRP_ALUM, I_COLLEC_DATA, rc)
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


  !! Calculates the emissions for CARMA aerosol particles. For alumina, the emissions
  !! are the discrete events that fire in this timestep (see alumina_events_mod). There
  !! is no surface flux.
  !!
  !! @author  Chuck Bardeen
  !! @version Jan-2011
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use ppgrid,             only: pcols, pver
    use physics_types,      only: physics_state
    use camsrfexch,         only: cam_in_t
    use alumina_events_mod, only: alumina_events_tendency

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

    ! Default return code.
    rc = RC_OK

    ncol = state%ncol

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8

    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.
    tendency(:, :) = 0.0_r8

    ! Add the events that fire in this timestep.
    if (ielem == I_ELEM_ALUM) then
      call alumina_events_tendency(state, ibin, dt, tendency)
    end if

    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use constituents,       only: pcnst
    use alumina_psd_mod,    only: alumina_psd_init
    use alumina_events_mod, only: alumina_events_init

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    real(kind=f)                       :: rlow(NBIN)            ! bin lower edge radius (cm)
    real(kind=f)                       :: rup(NBIN)             ! bin upper edge radius (cm)
    real(r8)                           :: rlow_r8(NBIN)         ! bin lower edge radius (cm)
    real(r8)                           :: rup_r8(NBIN)          ! bin upper edge radius (cm)

    integer                            :: LUNOPRT               ! logical unit number for output
    logical                            :: do_print              ! do print output?

    ! Default return code.
    rc = RC_OK

    ! Add initialization here.
    call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")

    ! Initialize the size distribution and the emission events.
    if (carma_do_emission) then

      ! The bin edges are set up by CARMA_Initialize, which is called before this.
      call CARMAGROUP_Get(carma, I_GRP_ALUM, rc, rlow=rlow, rup=rup)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMAGROUP_Get failed.")

      rlow_r8(:) = real(rlow(:), r8)
      rup_r8(:)  = real(rup(:), r8)

      call alumina_psd_init(rlow_r8, rup_r8, rc)
      if (rc < 0) call endrun("CARMA_InitializeModel: alumina_psd_init failed.")

      call alumina_events_init(rc)
      if (rc < 0) call endrun("CARMA_InitializeModel: alumina_events_init failed.")

      if (masterproc) then
        if (do_print) write(LUNOPRT,*) 'CARMA_InitializeModel: Done with alumina emission setup.'
      end if
    endif

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
