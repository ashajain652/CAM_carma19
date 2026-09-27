!! This module provides the particle size distribution (PSD) used to split the mass
!! of each alumina emission event across the CARMA bins.
!!
!! The PSD is either read from a netCDF file (carma_emis_psd_type = 'file') or built
!! from a lognormal (carma_emis_psd_type = 'lognormal'). A PSD file has:
!!
!!   dimensions : psd_alt, psd_bin, psd_bin_edge (= psd_bin + 1)
!!   psd_alt(psd_alt)             altitude above ground (km), increasing
!!   r_edge(psd_bin_edge)         bin edge radius (nm), increasing
!!   mass_frac(psd_alt, psd_bin)  mass fraction in each input bin, each row sums to 1
!!
!! At initialization, each altitude row of the PSD is remapped onto the CARMA bin
!! edges (rlow, rup) by overlap in ln r, assuming that the mass is uniform in ln r
!! inside each input bin. Mass that falls outside of the CARMA bin range is reported
!! as clipped, and each row is then renormalized to 1.
!!
!! This module does not depend upon the CAM physics state. It only needs the CARMA
!! bin edges, the namelist settings and MPI to broadcast the file contents.
!!
!! @version Sep-2026
module alumina_psd_mod

  use shr_kind_mod,          only: r8 => shr_kind_r8
  use spmd_utils,            only: masterproc
  use cam_abortutils,        only: endrun
  use cam_logfile,           only: iulog
  use carma_model_flags_mod, only: carma_emis_psd_type, carma_emis_psd_file, &
                                   carma_emis_psd_allow_clip, carma_emis_rmode, carma_emis_sigma

#if ( defined SPMD )
  use mpishorthand
#endif

  implicit none

  private

  ! Declare the public methods.
  public alumina_psd_init
  public alumina_psd_binfrac

  ! Declare the public data.
  integer, public, protected            :: psd_nbin = 0     !! number of CARMA bins in the remapped PSD

  ! Maximum fraction of the PSD mass that may fall outside of the CARMA bins.
  real(r8), parameter                   :: MAX_CLIP = 0.01_r8

  ! The PSD remapped onto the CARMA bins.
  integer                               :: psd_nalt = 0     ! number of PSD altitudes
  real(r8), allocatable, dimension(:)   :: psd_alt          ! PSD altitudes (km)
  real(r8), allocatable, dimension(:,:) :: psd_frac         ! mass fraction in each CARMA bin (psd_nbin, psd_nalt)

contains


  !! Initialize the PSD: read the file or build the lognormal, remap it onto the
  !! CARMA bins, check the clipped fraction and renormalize.
  !!
  !! @version Sep-2026
  subroutine alumina_psd_init(rlow, rup, rc)
    real(r8), intent(in)                :: rlow(:)        !! CARMA bin lower edge radius (cm)
    real(r8), intent(in)                :: rup(:)         !! CARMA bin upper edge radius (cm)
    integer, intent(out)                :: rc             !! return code, negative indicates failure

    integer                             :: nbin_in        ! number of input bins
    real(r8), allocatable               :: r_edge(:)      ! input bin edges (nm)
    real(r8), allocatable               :: mass_frac(:,:) ! input mass fractions (nbin_in, psd_nalt)
    real(r8), allocatable               :: clip(:)        ! clipped fraction at each altitude
    real(r8)                            :: lnlo(size(rlow)) ! ln of CARMA bin lower edge (nm)
    real(r8)                            :: lnhi(size(rlow)) ! ln of CARMA bin upper edge (nm)
    real(r8)                            :: rm             ! lognormal mass median radius (nm)
    real(r8)                            :: lnsig          ! ln of the geometric standard deviation
    real(r8)                            :: rowsum         ! sum of a row
    real(r8)                            :: w              ! overlap weight
    integer                             :: ia, i, j

    ! Default return code.
    rc = 0

    psd_nbin = size(rlow)
    if (size(rup) /= psd_nbin) call endrun('alumina_psd_init: rlow and rup have different sizes.')

    ! CARMA bin edges are in cm, while the PSD is in nm (1 nm = 1e-7 cm).
    lnlo(:) = log(rlow(:) * 1.e7_r8)
    lnhi(:) = log(rup(:)  * 1.e7_r8)

    if (trim(carma_emis_psd_type) == 'file') then

      call alumina_psd_read(nbin_in, r_edge, mass_frac)

      allocate(psd_frac(psd_nbin, psd_nalt))
      psd_frac(:,:) = 0._r8

      ! Normalize the input rows. The converter already does this, so only a
      ! round off correction is expected here.
      do ia = 1, psd_nalt
        rowsum = sum(mass_frac(:, ia))
        if (rowsum <= 0._r8) then
          write(iulog,*) 'alumina_psd_init: ERROR - mass_frac row sums to ', rowsum, ' at ', psd_alt(ia), ' km'
          call endrun('alumina_psd_init: PSD row does not have a positive sum.')
        end if
        if (masterproc .and. (abs(rowsum - 1._r8) > 1.e-6_r8)) then
          write(iulog,*) 'alumina_psd_init: WARNING - mass_frac row at ', psd_alt(ia), &
                         ' km sums to ', rowsum, ', renormalizing.'
        end if
        mass_frac(:, ia) = mass_frac(:, ia) / rowsum
      end do

      ! Remap by overlap in ln r. For input bin j and CARMA bin i:
      !
      !   W(i,j) = overlap_lnr(i,j) / dlnr(j)
      !   psd_frac(i,a) = sum_j W(i,j) * mass_frac(j,a)
      do j = 1, nbin_in
        do i = 1, psd_nbin
          w = min(lnhi(i), log(r_edge(j+1))) - max(lnlo(i), log(r_edge(j)))
          if (w > 0._r8) then
            w = w / (log(r_edge(j+1)) - log(r_edge(j)))
            psd_frac(i, :) = psd_frac(i, :) + w * mass_frac(j, :)
          end if
        end do
      end do

      deallocate(r_edge)
      deallocate(mass_frac)

    else if (trim(carma_emis_psd_type) == 'lognormal') then

      if (carma_emis_rmode <= 0._r8) call endrun('alumina_psd_init: carma_emis_rmode must be > 0.')
      if (carma_emis_sigma <= 1._r8) call endrun('alumina_psd_init: carma_emis_sigma must be > 1.')

      ! A single row that is used at all altitudes.
      psd_nalt = 1
      allocate(psd_alt(psd_nalt))
      allocate(psd_frac(psd_nbin, psd_nalt))
      psd_alt(1) = 0._r8

      ! carma_emis_rmode is the number median radius. The mass distribution of a
      ! lognormal number distribution is lognormal with the same sigma and a
      ! median radius of rmode * exp(3 ln^2 sigma). The mass in each CARMA bin is
      ! integrated exactly with erf, so the mass outside the bins shows up as
      ! clipped below.
      lnsig = log(carma_emis_sigma)
      rm    = carma_emis_rmode * exp(3._r8 * lnsig**2)

      do i = 1, psd_nbin
        psd_frac(i, 1) = 0.5_r8 * (erf((lnhi(i) - log(rm)) / (sqrt(2._r8) * lnsig)) - &
                                   erf((lnlo(i) - log(rm)) / (sqrt(2._r8) * lnsig)))
      end do

    else
      call endrun('alumina_psd_init: unknown carma_emis_psd_type ' // trim(carma_emis_psd_type) // &
                  ', expected file or lognormal.')
    end if

    ! Check the clipped fraction and renormalize.
    allocate(clip(psd_nalt))

    do ia = 1, psd_nalt
      rowsum   = sum(psd_frac(:, ia))
      clip(ia) = 1._r8 - rowsum

      if (masterproc) then
        write(iulog,'(a,f10.3,a,es13.5)') 'alumina_psd_init: altitude (km) = ', psd_alt(ia), &
          ', clipped fraction = ', clip(ia)
      end if

      if (rowsum <= 0._r8) then
        call endrun('alumina_psd_init: none of the PSD mass falls within the CARMA bins.')
      end if

      if (clip(ia) > MAX_CLIP) then
        if (carma_emis_psd_allow_clip) then
          if (masterproc) write(iulog,*) 'alumina_psd_init: WARNING - more than 1% of the PSD mass ', &
            'is outside of the CARMA bins at ', psd_alt(ia), ' km (allowed by carma_emis_psd_allow_clip).'
        else
          write(iulog,*) 'alumina_psd_init: ERROR - clipped fraction ', clip(ia), ' at ', psd_alt(ia), ' km'
          call endrun('alumina_psd_init: more than 1% of the PSD mass is outside of the CARMA bins. ' // &
                      'Change the CARMA bins or set carma_emis_psd_allow_clip.')
        end if
      end if

      psd_frac(:, ia) = psd_frac(:, ia) / rowsum
    end do

    ! Print the remapped table: one row per CARMA bin, one column per altitude.
    if (masterproc) call alumina_psd_print(lnlo, lnhi)

    deallocate(clip)

    return
  end subroutine alumina_psd_init


  !! Read the PSD file on masterproc and broadcast it to all tasks. Sets psd_nalt
  !! and psd_alt, and returns the input bin edges and mass fractions.
  !!
  !! NOTE: netCDF declares mass_frac(psd_alt, psd_bin) in C order. The Fortran
  !! interface reverses the dimension order, so the array is read here as
  !! mass_frac(psd_bin, psd_alt). This is checked against the variable's dimension
  !! ids.
  !!
  !! @version Sep-2026
  subroutine alumina_psd_read(nbin_in, r_edge, mass_frac)
    use ioFileMod,  only: getfil
    use wrap_nf,    only: wrap_open, wrap_close, wrap_inq_dimid, wrap_inq_dimlen, &
                          wrap_inq_varid, wrap_inq_varndims, wrap_inq_vardimid, &
                          wrap_get_var_realx, handle_error
    use netcdf,     only: nf90_get_var, NF90_NOERR, NF90_NOWRITE

    integer, intent(out)                  :: nbin_in        !! number of input bins
    real(r8), allocatable, intent(inout)  :: r_edge(:)      !! input bin edges (nm)
    real(r8), allocatable, intent(inout)  :: mass_frac(:,:) !! input mass fractions (nbin_in, psd_nalt)

    character(len=256)                    :: pfile          ! local PSD file name
    integer                               :: fid            ! file id
    integer                               :: alt_did        ! psd_alt dimension id
    integer                               :: bin_did        ! psd_bin dimension id
    integer                               :: edge_did       ! psd_bin_edge dimension id
    integer                               :: nedge          ! number of bin edges
    integer                               :: vid            ! variable id
    integer                               :: ndims          ! number of dimensions of a variable
    integer                               :: dimids(2)      ! dimension ids of a variable
    integer                               :: ret            ! netCDF return code
    integer                               :: ia, j

    if (masterproc) then
      call getfil(carma_emis_psd_file, pfile)
      write(iulog,*) 'alumina_psd_read: Reading particle size distribution from ', trim(pfile)

      call wrap_open(pfile, NF90_NOWRITE, fid)

      call wrap_inq_dimid(fid, 'psd_alt', alt_did)
      call wrap_inq_dimlen(fid, alt_did, psd_nalt)
      call wrap_inq_dimid(fid, 'psd_bin', bin_did)
      call wrap_inq_dimlen(fid, bin_did, nbin_in)
      call wrap_inq_dimid(fid, 'psd_bin_edge', edge_did)
      call wrap_inq_dimlen(fid, edge_did, nedge)

      if (nedge /= nbin_in + 1) call endrun('alumina_psd_read: psd_bin_edge must equal psd_bin + 1.')
      if ((psd_nalt < 1) .or. (nbin_in < 1)) call endrun('alumina_psd_read: empty PSD file.')
    end if

#if ( defined SPMD )
    call mpibcast(psd_nalt, 1, mpiint, 0, mpicom)
    call mpibcast(nbin_in,  1, mpiint, 0, mpicom)
#endif

    allocate(psd_alt(psd_nalt))
    allocate(r_edge(nbin_in+1))
    allocate(mass_frac(nbin_in, psd_nalt))

    if (masterproc) then
      call wrap_inq_varid(fid, 'psd_alt', vid)
      call wrap_get_var_realx(fid, vid, psd_alt)

      call wrap_inq_varid(fid, 'r_edge', vid)
      call wrap_get_var_realx(fid, vid, r_edge)

      ! The file declares mass_frac(psd_alt, psd_bin); in Fortran (reversed) order
      ! the first dimension must be psd_bin and the second psd_alt.
      call wrap_inq_varid(fid, 'mass_frac', vid)
      call wrap_inq_varndims(fid, vid, ndims)
      if (ndims /= 2) call endrun('alumina_psd_read: mass_frac must have 2 dimensions.')
      call wrap_inq_vardimid(fid, vid, dimids)
      if ((dimids(1) /= bin_did) .or. (dimids(2) /= alt_did)) then
        call endrun('alumina_psd_read: mass_frac must be declared as mass_frac(psd_alt, psd_bin) in the file.')
      end if

      ret = nf90_get_var(fid, vid, mass_frac)
      if (ret /= NF90_NOERR) then
        write(iulog,*) 'alumina_psd_read: error reading mass_frac'
        call handle_error(ret)
      end if

      call wrap_close(fid)

      ! Sanity checks.
      do ia = 2, psd_nalt
        if (psd_alt(ia) <= psd_alt(ia-1)) call endrun('alumina_psd_read: psd_alt must be strictly increasing.')
      end do
      if (r_edge(1) <= 0._r8) call endrun('alumina_psd_read: r_edge must be positive.')
      do j = 2, nbin_in + 1
        if (r_edge(j) <= r_edge(j-1)) call endrun('alumina_psd_read: r_edge must be strictly increasing.')
      end do
      if (any(mass_frac < 0._r8)) call endrun('alumina_psd_read: mass_frac must be >= 0.')

      write(iulog,*) 'alumina_psd_read: psd_nalt = ', psd_nalt, ', psd_bin = ', nbin_in
      write(iulog,*) 'alumina_psd_read: r_edge (nm) from ', r_edge(1), ' to ', r_edge(nbin_in+1)
      write(iulog,*) 'alumina_psd_read: psd_alt (km) from ', psd_alt(1), ' to ', psd_alt(psd_nalt)
    end if

#if ( defined SPMD )
    call mpibcast(psd_alt,   psd_nalt,           mpir8, 0, mpicom)
    call mpibcast(r_edge,    nbin_in+1,          mpir8, 0, mpicom)
    call mpibcast(mass_frac, nbin_in*psd_nalt,   mpir8, 0, mpicom)
#endif

    return
  end subroutine alumina_psd_read


  !! Print the remapped PSD table to the log.
  !!
  !! @version Sep-2026
  subroutine alumina_psd_print(lnlo, lnhi)
    real(r8), intent(in)                :: lnlo(:)      !! ln of CARMA bin lower edge (nm)
    real(r8), intent(in)                :: lnhi(:)      !! ln of CARMA bin upper edge (nm)

    character(len=64)                   :: fmt          ! format string
    integer                             :: i

    write(iulog,*) ''
    write(iulog,*) 'alumina_psd_init: PSD mass fraction (psdfrac) in each CARMA bin, one column per altitude (km)'
    write(fmt, '(a,i0,a)') '(a4,2a12,', psd_nalt, 'f14.3)'
    write(iulog, fmt) 'bin', 'rlow (nm)', 'rup (nm)', psd_alt(:)
    write(fmt, '(a,i0,a)') '(i4,2es12.4,', psd_nalt, 'es14.6)'
    do i = 1, psd_nbin
      write(iulog, fmt) i, exp(lnlo(i)), exp(lnhi(i)), psd_frac(i, :)
    end do
    write(iulog,*) ''

    return
  end subroutine alumina_psd_print


  !! Return the fraction of the emitted mass that goes into each CARMA bin at the
  !! altitude z_km. The PSD is interpolated linearly between the PSD altitudes and
  !! clamped outside of their range.
  !!
  !! @version Sep-2026
  subroutine alumina_psd_binfrac(z_km, binfrac)
    real(r8), intent(in)                :: z_km          !! altitude above ground (km)
    real(r8), intent(out)               :: binfrac(:)    !! mass fraction in each CARMA bin

    integer                             :: ia
    real(r8)                            :: w             ! interpolation weight

    if (.not. allocated(psd_frac)) call endrun('alumina_psd_binfrac: PSD is not initialized.')
    if (size(binfrac) /= psd_nbin) call endrun('alumina_psd_binfrac: binfrac has the wrong size.')

    if ((psd_nalt == 1) .or. (z_km <= psd_alt(1))) then
      binfrac(:) = psd_frac(:, 1)
    else if (z_km >= psd_alt(psd_nalt)) then
      binfrac(:) = psd_frac(:, psd_nalt)
    else
      do ia = 1, psd_nalt - 1
        if (z_km < psd_alt(ia+1)) exit
      end do

      ! psd_alt(ia) <= z_km < psd_alt(ia+1)
      w = (z_km - psd_alt(ia)) / (psd_alt(ia+1) - psd_alt(ia))
      binfrac(:) = (1._r8 - w) * psd_frac(:, ia) + w * psd_frac(:, ia+1)
    end if

    return
  end subroutine alumina_psd_binfrac

end module alumina_psd_mod
