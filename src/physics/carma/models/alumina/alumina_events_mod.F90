!! This module handles the discrete emission events for the alumina model. Each
!! row of the event file is a point release of a mass of Al2O3 (kg) at a location
!! and an altitude above ground, at a date and time in the model calendar. All of a
!! row's mass is injected within one model timestep (a dirac pulse) into the physics
!! column nearest to the row and the layer that contains the row's altitude. The mass
!! is split across the CARMA bins using the PSD at the row's altitude (see
!! alumina_psd_mod).
!!
!! The event file (netCDF) has:
!!
!!   dimension : row
!!   event_id(row)  int    event identifier (shared by the rows of one event)
!!   date(row)      int    date (YYYYMMDD)
!!   datesec(row)   int    seconds of the day
!!   lat(row)       real   latitude (degrees north)
!!   lon(row)       real   longitude (degrees east, 0-360)
!!   altitude(row)  real   altitude above ground (km)
!!   mass(row)      real   emitted mass (kg)
!!
!! @version Sep-2026
module alumina_events_mod

  use shr_kind_mod,          only: r8 => shr_kind_r8
  use spmd_utils,            only: masterproc
  use cam_abortutils,        only: endrun
  use cam_logfile,           only: iulog
  use carma_model_flags_mod, only: carma_emis_file, carma_emis_scale
  use alumina_psd_mod,       only: psd_nbin, alumina_psd_binfrac

#if ( defined SPMD )
  use mpishorthand
#endif

  implicit none

  private

  ! Declare the public methods.
  public alumina_events_init
  public alumina_events_tendency

  ! Event table, identical on all tasks.
  integer                               :: evt_nrow = 0   ! number of rows
  integer, allocatable, dimension(:)    :: evt_id         ! event id
  integer, allocatable, dimension(:)    :: evt_date       ! date (YYYYMMDD)
  integer, allocatable, dimension(:)    :: evt_sec        ! seconds of the day
  real(r8), allocatable, dimension(:)   :: evt_lat        ! latitude (degrees)
  real(r8), allocatable, dimension(:)   :: evt_lon        ! longitude (degrees, 0-360)
  real(r8), allocatable, dimension(:)   :: evt_alt        ! altitude above ground (km)
  real(r8), allocatable, dimension(:)   :: evt_mass       ! mass (kg)
  integer, allocatable, dimension(:)    :: evt_gcol       ! global column that receives the row
  real(r8), allocatable, dimension(:,:) :: evt_binfrac    ! mass fraction in each CARMA bin (evt_nrow, psd_nbin)

  ! Location of the receiving column on this task, or -1 if another task owns it.
  integer, allocatable, dimension(:)    :: evt_lchnk      ! local chunk index
  integer, allocatable, dimension(:)    :: evt_icol       ! column index in the chunk

contains


  !! Read the event file, map each row to the nearest physics column and compute
  !! the size split of each row.
  !!
  !! NOTE: alumina_psd_init must be called first.
  !!
  !! @version Sep-2026
  subroutine alumina_events_init(rc)
    use ioFileMod,     only: getfil
    use wrap_nf,       only: wrap_open, wrap_close, wrap_inq_dimid, wrap_inq_dimlen, &
                             wrap_inq_varid, wrap_get_var_realx, wrap_get_var_int
    use netcdf,        only: NF90_NOWRITE
    use ppgrid,        only: pcols, begchunk, endchunk
    use phys_grid,     only: get_ncols_p, get_rlat_all_p, get_rlon_all_p, get_gcol_p
    use physconst,     only: pi, rearth
    use time_manager,  only: get_prev_date

    integer, intent(out)                :: rc             !! return code, negative indicates failure

    character(len=256)                  :: efile          ! local event file name
    integer                             :: fid            ! file id
    integer                             :: row_did        ! row dimension id
    integer                             :: vid            ! variable id
    integer                             :: irow           ! row index
    integer                             :: lchnk          ! chunk index
    integer                             :: ncol           ! number of columns in the chunk
    integer                             :: icol           ! column index
    integer                             :: nlate          ! number of rows before the current time
    integer                             :: yr, mon, day, tod, ymd
    real(r8)                            :: rlats(pcols)   ! column latitudes (radians)
    real(r8)                            :: rlons(pcols)   ! column longitudes (radians)
    real(r8)                            :: elat, elon     ! row latitude and longitude (radians)
    real(r8)                            :: dist           ! great circle distance (radians)
    real(r8)                            :: deg2rad        ! degrees to radians
    real(r8), allocatable               :: dmin(:,:)      ! local (distance, gcol) pairs (2, evt_nrow)
    real(r8), allocatable               :: gmin(:,:)      ! global (distance, gcol) pairs (2, evt_nrow)
    real(r8), allocatable               :: lclatlon(:,:)  ! local receiving column lat/lon (2, evt_nrow)
    real(r8), allocatable               :: glatlon(:,:)   ! global receiving column lat/lon (2, evt_nrow)
    integer, allocatable                :: my_lchnk(:)    ! local candidate chunk
    integer, allocatable                :: my_icol(:)     ! local candidate column
#if ( defined SPMD )
    integer                             :: ier            ! MPI error code
#endif

    ! Default return code.
    rc = 0

    deg2rad = pi / 180._r8

    if (psd_nbin < 1) call endrun('alumina_events_init: alumina_psd_init must be called first.')

    ! Read the event file.
    if (masterproc) then
      call getfil(carma_emis_file, efile)
      write(iulog,*) 'alumina_events_init: Reading emission events from ', trim(efile)

      call wrap_open(efile, NF90_NOWRITE, fid)
      call wrap_inq_dimid(fid, 'row', row_did)
      call wrap_inq_dimlen(fid, row_did, evt_nrow)
    end if

#if ( defined SPMD )
    call mpibcast(evt_nrow, 1, mpiint, 0, mpicom)
#endif

    if (evt_nrow < 1) then
      if (masterproc) then
        call wrap_close(fid)
        write(iulog,*) 'alumina_events_init: WARNING - the event file has no rows, no emissions.'
      end if
      return
    end if

    allocate(evt_id(evt_nrow), evt_date(evt_nrow), evt_sec(evt_nrow))
    allocate(evt_lat(evt_nrow), evt_lon(evt_nrow), evt_alt(evt_nrow), evt_mass(evt_nrow))
    allocate(evt_gcol(evt_nrow), evt_lchnk(evt_nrow), evt_icol(evt_nrow))
    allocate(evt_binfrac(evt_nrow, psd_nbin))

    if (masterproc) then
      call wrap_inq_varid(fid, 'event_id', vid)
      call wrap_get_var_int(fid, vid, evt_id)
      call wrap_inq_varid(fid, 'date', vid)
      call wrap_get_var_int(fid, vid, evt_date)
      call wrap_inq_varid(fid, 'datesec', vid)
      call wrap_get_var_int(fid, vid, evt_sec)
      call wrap_inq_varid(fid, 'lat', vid)
      call wrap_get_var_realx(fid, vid, evt_lat)
      call wrap_inq_varid(fid, 'lon', vid)
      call wrap_get_var_realx(fid, vid, evt_lon)
      call wrap_inq_varid(fid, 'altitude', vid)
      call wrap_get_var_realx(fid, vid, evt_alt)
      call wrap_inq_varid(fid, 'mass', vid)
      call wrap_get_var_realx(fid, vid, evt_mass)

      call wrap_close(fid)

      ! Sanity checks.
      do irow = 1, evt_nrow
        if ((evt_lat(irow) < -90._r8) .or. (evt_lat(irow) > 90._r8) .or. &
            (evt_sec(irow) < 0) .or. (evt_sec(irow) >= 86400) .or. &
            (evt_mass(irow) < 0._r8) .or. (evt_alt(irow) < 0._r8)) then
          write(iulog,*) 'alumina_events_init: ERROR - bad row ', irow, ' event ', evt_id(irow), &
            evt_date(irow), evt_sec(irow), evt_lat(irow), evt_lon(irow), evt_alt(irow), evt_mass(irow)
          call endrun('alumina_events_init: bad row in the event file.')
        end if
      end do

      write(iulog,*) 'alumina_events_init: rows = ', evt_nrow, ', total mass (kg) = ', sum(evt_mass), &
        ', carma_emis_scale = ', carma_emis_scale
    end if

#if ( defined SPMD )
    call mpibcast(evt_id,   evt_nrow, mpiint, 0, mpicom)
    call mpibcast(evt_date, evt_nrow, mpiint, 0, mpicom)
    call mpibcast(evt_sec,  evt_nrow, mpiint, 0, mpicom)
    call mpibcast(evt_lat,  evt_nrow, mpir8,  0, mpicom)
    call mpibcast(evt_lon,  evt_nrow, mpir8,  0, mpicom)
    call mpibcast(evt_alt,  evt_nrow, mpir8,  0, mpicom)
    call mpibcast(evt_mass, evt_nrow, mpir8,  0, mpicom)
#endif

    ! Find the nearest physics column on this task for each row, by great circle
    ! distance (haversine). The pair (distance, gcol) is kept for a MINLOC reduction.
    allocate(dmin(2, evt_nrow), gmin(2, evt_nrow))
    allocate(lclatlon(2, evt_nrow), glatlon(2, evt_nrow))
    allocate(my_lchnk(evt_nrow), my_icol(evt_nrow))

    dmin(1, :) = huge(1._r8)
    dmin(2, :) = huge(1._r8)
    my_lchnk(:) = -1
    my_icol(:)  = -1

    do lchnk = begchunk, endchunk
      ncol = get_ncols_p(lchnk)
      call get_rlat_all_p(lchnk, pcols, rlats)
      call get_rlon_all_p(lchnk, pcols, rlons)

      do irow = 1, evt_nrow
        elat = evt_lat(irow) * deg2rad
        elon = evt_lon(irow) * deg2rad

        do icol = 1, ncol
          dist = 2._r8 * asin(min(1._r8, sqrt(sin(0.5_r8 * (rlats(icol) - elat))**2 + &
                 cos(elat) * cos(rlats(icol)) * sin(0.5_r8 * (rlons(icol) - elon))**2)))

          ! Ties are broken by the lowest global column index, the same as MINLOC.
          if ((dist < dmin(1, irow)) .or. &
              ((dist == dmin(1, irow)) .and. (real(get_gcol_p(lchnk, icol), r8) < dmin(2, irow)))) then
            dmin(1, irow) = dist
            dmin(2, irow) = real(get_gcol_p(lchnk, icol), r8)
            my_lchnk(irow) = lchnk
            my_icol(irow)  = icol
          end if
        end do
      end do
    end do

    ! Pick the owning column over all tasks. MPI_MINLOC on (value, index) pairs
    ! returns the smallest distance, and the smallest gcol if distances are equal.
#if ( defined SPMD )
    call mpi_allreduce(dmin, gmin, evt_nrow, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mpicom, ier)
    if (ier /= MPI_SUCCESS) call endrun('alumina_events_init: mpi_allreduce (MINLOC) failed.')
#else
    gmin(:, :) = dmin(:, :)
#endif

    ! Keep the rows whose column is on this task.
    lclatlon(:, :) = 0._r8
    do irow = 1, evt_nrow
      evt_gcol(irow)  = nint(gmin(2, irow))
      evt_lchnk(irow) = -1
      evt_icol(irow)  = -1

      if ((my_lchnk(irow) > 0) .and. (nint(dmin(2, irow)) == evt_gcol(irow))) then
        evt_lchnk(irow) = my_lchnk(irow)
        evt_icol(irow)  = my_icol(irow)

        call get_rlat_all_p(my_lchnk(irow), pcols, rlats)
        call get_rlon_all_p(my_lchnk(irow), pcols, rlons)
        lclatlon(1, irow) = rlats(my_icol(irow)) / deg2rad
        lclatlon(2, irow) = rlons(my_icol(irow)) / deg2rad
      end if
    end do

    ! Gather the location of the receiving columns for the log. Only the owner
    ! contributes a nonzero value.
#if ( defined SPMD )
    call mpi_allreduce(lclatlon, glatlon, 2*evt_nrow, MPI_REAL8, MPI_SUM, mpicom, ier)
    if (ier /= MPI_SUCCESS) call endrun('alumina_events_init: mpi_allreduce (SUM) failed.')
#else
    glatlon(:, :) = lclatlon(:, :)
#endif

    ! Precompute the size split of each row at the row's own altitude.
    do irow = 1, evt_nrow
      call alumina_psd_binfrac(evt_alt(irow), evt_binfrac(irow, :))
    end do

    ! Log the mapping.
    if (masterproc) then
      write(iulog,*) ''
      write(iulog,*) 'alumina_events_init: row, event id, date, sec, altitude (km), mass (kg), ', &
        'lat, lon, gcol, column lat, column lon, distance (km)'
      do irow = 1, evt_nrow
        write(iulog,'(a,i6,i8,i10,i7,f9.3,es14.6,2f10.4,i8,2f10.4,f10.2)') ' alumina_events_init: ', &
          irow, evt_id(irow), evt_date(irow), evt_sec(irow), evt_alt(irow), evt_mass(irow), &
          evt_lat(irow), evt_lon(irow), evt_gcol(irow), glatlon(1, irow), glatlon(2, irow), &
          gmin(1, irow) * rearth * 1.e-3_r8
      end do
      write(iulog,*) ''

      ! Rows before the start of the current timestep will never fire in this run
      ! segment. NOTE: The CESM time manager defines the previous time as the current
      ! time minus one timestep, even on the first step.
      call get_prev_date(yr, mon, day, tod)
      ymd = yr*10000 + mon*100 + day
      nlate = 0
      do irow = 1, evt_nrow
        if (time_lt(evt_date(irow), evt_sec(irow), ymd, tod)) nlate = nlate + 1
      end do
      if (nlate > 0) then
        write(iulog,*) 'alumina_events_init: WARNING - ', nlate, ' rows are before the model time ', &
          ymd, tod, ' and will not fire in this run segment.'
      end if
    end if

    deallocate(dmin, gmin, lclatlon, glatlon, my_lchnk, my_icol)

    return
  end subroutine alumina_events_init


  !! Add the tendency (kg/kg/s) for bin ibin from all rows that fire in this timestep
  !! and whose receiving column is in this chunk. A row fires when
  !!
  !!   prev_date (start of step) <= row time < curr_date (end of step)
  !!
  !! so each row fires exactly once. Dates are compared as (YYYYMMDD, sec) integer
  !! pairs, which are ordered the same as time in any calendar.
  !!
  !! NOTE: tendency is not reset here; the caller must initialize it.
  !!
  !! @version Sep-2026
  subroutine alumina_events_tendency(state, ibin, dt, tendency)
    use physics_types, only: physics_state
    use ppgrid,        only: pcols, pver, pverp
    use phys_grid,     only: get_area_p
    use physconst,     only: pi, rearth, gravit
    use time_manager,  only: get_prev_date, get_curr_date, get_nstep

    type(physics_state), intent(in)     :: state                 !! physics state
    integer, intent(in)                 :: ibin                  !! bin index
    real(r8), intent(in)                :: dt                    !! time step (s)
    real(r8), intent(inout)             :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)

    integer                             :: irow                  ! row index
    integer                             :: icol                  ! column index
    integer                             :: k                     ! layer index
    integer                             :: kfound                ! layer that contains the altitude
    integer                             :: lchnk                 ! chunk index
    integer                             :: yr, mon, day
    integer                             :: prev_ymd, prev_tod    ! start of the timestep
    integer                             :: curr_ymd, curr_tod    ! end of the timestep
    real(r8)                            :: alt_m                 ! altitude above ground (m)
    real(r8)                            :: area                  ! column area (m2)
    real(r8)                            :: airmass               ! mass of air in the layer (kg)
    real(r8)                            :: mass                  ! mass emitted into this bin (kg)

    if (evt_nrow < 1) return

    lchnk = state%lchnk

    ! Nothing to do unless a row is in this chunk.
    if (.not. any(evt_lchnk(:) == lchnk)) return

    call get_prev_date(yr, mon, day, prev_tod)
    prev_ymd = yr*10000 + mon*100 + day
    call get_curr_date(yr, mon, day, curr_tod)
    curr_ymd = yr*10000 + mon*100 + day

    do irow = 1, evt_nrow
      if (evt_lchnk(irow) /= lchnk) cycle

      ! prev <= t_event < curr
      if (time_lt(evt_date(irow), evt_sec(irow), prev_ymd, prev_tod)) cycle
      if (.not. time_lt(evt_date(irow), evt_sec(irow), curr_ymd, curr_tod)) cycle

      icol  = evt_icol(irow)
      alt_m = evt_alt(irow) * 1000._r8

      ! Find the layer with zi(k+1) <= alt < zi(k). zi is the height above the
      ! surface at the interfaces (m), with k = 1 at the model top.
      kfound = -1
      do k = 1, pver
        if ((state%zi(icol, k+1) <= alt_m) .and. (alt_m < state%zi(icol, k))) then
          kfound = k
          exit
        end if
      end do

      if (kfound < 0) then
        write(iulog,*) 'alumina_events_tendency: ERROR - row ', irow, ' event ', evt_id(irow), &
          ' altitude (m) ', alt_m, ' is outside of the column, zi(top) = ', state%zi(icol, 1), &
          ', zi(surface) = ', state%zi(icol, pverp)
        call endrun('alumina_events_tendency: emission altitude is outside of the model column.')
      end if

      ! Convert the mass to a tendency on the (wet) mass mixing ratio.
      area    = get_area_p(lchnk, icol) * rearth**2
      airmass = area * state%pdel(icol, kfound) / gravit
      mass    = evt_mass(irow) * carma_emis_scale * evt_binfrac(irow, ibin)

      tendency(icol, kfound) = tendency(icol, kfound) + mass / (dt * airmass)

      if (ibin == 1) then
        write(iulog,'(a,i8,a,i6,a,i8,a,i10,i7,a,f9.3,a,2f10.4,a,i8,a,i4,a,es14.6)') &
          'alumina_events_tendency: nstep ', get_nstep(), ' row ', irow, ' event ', evt_id(irow), &
          ' date ', evt_date(irow), evt_sec(irow), ' alt (km) ', evt_alt(irow), &
          ' col lat/lon ', state%lat(icol) * 180._r8 / pi, state%lon(icol) * 180._r8 / pi, &
          ' gcol ', evt_gcol(irow), ' k ', kfound, ' mass (kg) ', evt_mass(irow) * carma_emis_scale
      end if
    end do

    return
  end subroutine alumina_events_tendency


  !! Returns true if (ymd1, tod1) is earlier than (ymd2, tod2).
  !!
  !! @version Sep-2026
  logical function time_lt(ymd1, tod1, ymd2, tod2)
    integer, intent(in)                 :: ymd1   !! date 1 (YYYYMMDD)
    integer, intent(in)                 :: tod1   !! seconds of day 1
    integer, intent(in)                 :: ymd2   !! date 2 (YYYYMMDD)
    integer, intent(in)                 :: tod2   !! seconds of day 2

    time_lt = (ymd1 < ymd2) .or. ((ymd1 == ymd2) .and. (tod1 < tod2))

    return
  end function time_lt

end module alumina_events_mod
