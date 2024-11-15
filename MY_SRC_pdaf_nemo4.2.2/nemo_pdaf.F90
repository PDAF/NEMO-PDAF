!> Include NEMO variables and initialize NEMO grid information
!!
!! This module includes variables from NEMO. For PDAF in its
!! online coupling it is the single point which directly links
!! to NEMO. All other routines include from this module.
!! For the offline case the variables are declared here and
!! separately initialized.
!!
!! Next to including information from NEMO, the routine
!! set_nemo_grid initializes arrays holding the grid information
!! for use with state vectors in the PDAF user code for NEMO.
!!
module nemo_pdaf

  use mod_kind_pdaf

#ifndef PDAF_OFFLINE
  ! Include variables from NEMO
  ! User routines should only include from `nemo_pdaf` not from NEMO modules
  ! The only exception is `asminc_pdaf` which also directly includes from NEMO
  use par_oce, &
       only: jpi, jpj, jpk, jpiglo, jpjglo, &
       Nis0, Nie0, Njs0, Nje0, nn_hls, &
       jp_tem, jp_sal
  use dom_oce, &
       only: glamt, gphit, nimpp, njmpp, gdept_1d, &
       ndastp, l_1st_euler, tmask, umask, vmask
  use oce, &
       only: ssh, ts, uu, vv
  use in_out_manager, &
       only: nitend, nit000, lwp, numout
  use lbclnk, &
       only: lbc_lnk
  use diaobs, &
       only: calc_date
#if defined key_top
  use trc, &
       only: tr
  use par_trc, &
       only: jptra
  use trcnam, &
       only: sn_tracer
#endif
#endif

  implicit none

  integer :: Nbb, Nnn, Nrhs

  ! *** NEmo grid dimensions - generic names following naming scheme of NEMO 4.0
  integer :: nldi, nldj                 ! first inner index in i/j direction of sub-domain
  integer :: nlei, nlej                 ! last inner index in i/j direction of sub-domain

#if defined PDAF_OFFLINE
  ! *** NEMO model variables - they used in the offline mode
  integer :: jpiglo, jpjglo, jpk        ! Global NEMO grid dimensions
  integer :: nimpp, njmpp               ! start i,j of subdomain including halo

  real(pwp), allocatable :: glamt(:,:)       ! Longitudes
  real(pwp), allocatable :: gphit(:,:)       ! Latitudes
  real(pwp), allocatable :: gdept_1d(:)      ! Depths
  real(pwp), allocatable :: tmask(:,:,:)     ! Temperature mask array

  logical :: lwp = .true.
  integer :: numout = 6
  integer :: jptra = 0
  integer :: ndastp = 1

#if defined key_top
  type tracer
     character(len=lc) :: clsname
     character(len=lc) :: clunit
  end type tracer
  
  type(tracer), allocatable :: sn_tracer(:)
#endif
#endif

  ! *** Other grid variables
  real(pwp), allocatable :: tmp_4d(:,:,:,:)     ! 4D array used to represent full NEMO grid box
  real(4), allocatable :: stmp_4d(:,:,:,:)      ! 4D array used to represent full NEMO grid box
  real(pwp), allocatable :: lat1(:), lon1(:)    ! Vectors holding latitude and latitude

  integer :: dim_2d                        ! Dimension of 2d grid box
  integer :: nwet                          ! Number of surface wet grid points
  integer :: nwet3d                        ! Number of 3d wet grid points
  integer, allocatable :: wet_pts(:,:)     ! Index array for wet grid points
                           ! (1) latitude, (2) langitude, (3) number wet layers, (4) index in 2d grid box
  integer, allocatable :: idx_wet_2d(:,:)  ! Index array for wet_pts row index in 2d box
  integer, allocatable :: idx_nwet(:,:)    ! Index array for wet_pts row index in wet surface grid points
  integer, allocatable :: nlev_wet_2d(:,:) ! Number of wet layers for ij position in 2d box

  integer :: use_wet_state=0               ! 1: State vector contains full columns where surface grid point is wet
                                           ! 2: State vector only contains wet grid points
                                           ! other: State vector contains 2d/3d grid box

  integer :: i0, j0                         ! PE-local halo offsets
  integer :: ni_p, nj_p, nk_p               ! Size of decomposed grid
  integer :: istart, jstart                 ! Start indices for internal local domain
  integer :: dim_2d_p, dim_3d_p             ! Dimension of 2d/3d grid box of sub-domain
  integer :: type_limcoords = 0             ! How limiting domain coords are determined
       ! 0: from glamt/gphit; 1: from min/max of lat1_p and lon1_p
  real(pwp), allocatable :: lat1_p(:), lon1_p(:) ! Vectors holding latitude and latitude for decomposition
  real(pwp), allocatable :: lats(:,:), lons(:,:) ! Arrays for interior coordinates (no halo)
  integer :: sdim2d, sdim3d                 ! 2D/3D dimension of field in state vector

  ! *** File name and path to read grid information
  character(len=200)  :: path_dims         ! Path for NEMO file holding dimensions
  character(len=80)   :: file_dims         ! File name NEMO file holding dimensions

! Constants for coordinate calculations
  real(8), parameter  :: pi = 3.14159265358979323846_pwp   ! Pi
  real :: deg2rad = pi / 180.0_pwp      ! Conversion from degrees to radian


contains

  !> Initialize grid information for the DA
  !!
  !! This routine initializes grid information for the data assimilation
  !! In particular index information is initialize to map in between
  !! the model grid and the state vector
  !!
  subroutine set_nemo_grid(screen)

    use mod_kind_pdaf
    use parallel_pdaf, &
         only: mype_model, npes_model, task_id, comm_model, &
         MPI_INT, MPI_SUM, MPIerr
    use PDAFomi, &
         only: PDAFomi_set_domain_limits

    implicit none

! *** Argument ***
    integer, intent(in) :: screen         ! Control verbosity

! *** Local variables ***
    integer :: i, j, k                    ! Counters
    integer :: cnt, cnt_all, cnt_layers   ! Counters
    integer :: nwet_g, nwet3d_g           ! Global sums of wet grid point
    real(pwp) :: lim_coords(2,2)          ! Limiting coordinates of sub-domain


! *** Get index dimensions from step or stpmfl

    ! This is in a subroutine to avoid cyclic dependency of nemo_pdaf and assimilation_pdaf
    CALL get_nemo_indices(Nbb, Nnn, Nrhs)

! *** set dimension of 2d and 3d fields in state vector ***

    ! Store generic names (following NEMO 4.0 naming)
    nlei = Nie0
    nlej = Nje0
    nldi = Nis0
    nldj = Njs0

    ! Local dimensions
    ni_p = nlei - nldi + 1
    nj_p = nlej - nldj + 1
    nk_p = jpk

    ! Compute halo offset
    i0 = nldi - nn_hls
    j0 = nldj - nn_hls

    ! Size of 2d/3d boxes without halo
    dim_2d_p = ni_p * nj_p
    dim_3d_p = ni_p * nj_p * nk_p

    ! Start indices for sub-domain without halo
    istart = nimpp + nldi - 1
    jstart = njmpp + nldj - 1


#ifndef PDAF_OFFLINE
    ! Set coordinate vectors for rectangular grids
    allocate(lat1_p(nj_p), lon1_p(ni_p))

    lat1_p(:) = 0.0
    do j = 1, nj_p
       do i = 1, ni_p
          if (abs(gphit(i+i0, j+j0)) > 0.00001) then
             lat1_p(j) = gphit(i+i0, j+j0)
          endif
       enddo
    enddo

    lon1_p(:) = 0.0
    do j = 1, nj_p
       do i = 1, ni_p
          if (abs(glamt(i+i0, j+j0)) > 0.00001) then
             lon1_p(i) = glamt(i+i0, j+j0)
          endif
       enddo
    enddo
#endif

    ! Store interior coordinates
    allocate(lons(ni_p, nj_p), lats(ni_p, nj_p))
    do j = 1, nj_p
       do i = 1, ni_p
          lats(i,j) = gphit(i0+i, j0+j)
          lons(i,j) = glamt(i0+i, j0+j)
       end do
    end do

    ! Count number of surface points
    cnt = 0
    do k = 1, nk_p
       do j = 1, nj_p
          do i = 1, ni_p 
             cnt = cnt + 1
             if (tmask(i + i0, j + j0, k) == 1.0_pwp) then
                if (k==1) nwet = nwet + 1
                nwet3d = nwet3d + 1
             endif
          enddo
       enddo
    enddo

    ! Initialize index arrays
    ! - for mapping from nx*ny grid to vector of wet points
    ! - mask for wet points

    if (nwet > 0) then

       allocate(wet_pts(7, nwet))

       cnt = 0
       cnt_all = 0
       do j = 1, nj_p
          do i = 1, ni_p
             cnt_all = cnt_all + 1
             if (tmask(i + i0, j + j0, 1) == 1.0_pwp) then
                cnt = cnt + 1

                wet_pts(1,cnt) = i + nimpp + nldi - 2     ! Global longitude index
                wet_pts(2,cnt) = j + njmpp + nldj - 2     ! Global latitude index
                wet_pts(6,cnt) = i                ! Longitude index in subdomain
                wet_pts(7,cnt) = j                ! Latitidue index in subdomain

                ! Determine number of wet layers
                cnt_layers = 0
                do k = 1, nk_p
                   if (tmask(i + i0, j + j0, k) == 1.0_pwp) cnt_layers = cnt_layers + 1
                end do
                wet_pts(3,cnt) = cnt_layers
                wet_pts(4,cnt) = cnt_all
             end if
          end do
       end do

       ! row 5 stores wet_pts index for vertical column storage in state vector
       wet_pts(5,1) = 1
       do i = 2 , nwet
          wet_pts(5,i) = wet_pts(5,i-1) + wet_pts(3,i-1)
       end do

    else

       write (*, '(8x,a,i3)') 'WARNING: No valid local domains, PE=', mype_model
       nwet = 0

       allocate(wet_pts(3, 1))

    end if

  ! Initialize index arrays
  ! - for mapping from vector of wet points to 2d box
  ! - for mapping from vector model grid point coordinats to wet point index
  ! - for retrieving number of wet layers at a grid point coordinate
  ! these arrays also serve as mask arrays (0 indicates land point)

  allocate(idx_wet_2d(ni_p, nj_p))
  allocate(idx_nwet(ni_p, nj_p))
  allocate(nlev_wet_2d(ni_p, nj_p))

  idx_wet_2d = 0
  idx_nwet = 0
  nlev_wet_2d = 0
  do i = 1 , nwet
     idx_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(4,i)
     nlev_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(3,i)
  end do
  if (use_wet_state/=2) then
     do i = 1 , nwet
        idx_nwet(wet_pts(6,i), wet_pts(7,i)) = i
     end do
  else
     do i = 1 , nwet
        idx_nwet(wet_pts(6,i), wet_pts(7,i)) = wet_pts(5,i)
     end do
  end if


! ********************************************************
! *** Set dimension of 2D and 3D field in state vector ***
! ********************************************************

    if (use_wet_state==1) then
       ! State vector contains full columns when surface grid point is wet
       sdim3d = abs(nwet)*jpk
       sdim2d = abs(nwet)
    elseif (use_wet_state==2) then
       ! State vector only contains wet grid points
       sdim3d = abs(nwet3d)
       sdim2d = abs(nwet)
    else
       ! State vector contains 2d/3d grid box
       sdim3d = dim_3d_p
       sdim2d = dim_2d_p
    end if

! *** Screen output ***

    if (npes_model==1) then
       write(*,'(a,5x,a,3x,i11)') 'NEMO-PDAF', 'Number of wet surface points', nwet
       write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', 'Number of 3D wet points', nwet3d
       write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', '2D wet points * nlayers', nwet*nk_p
    else 
       if (screen==1 .or. screen==2) then

          if (task_id==1) then
             ! Get global sums
             call MPI_Reduce (nwet, nwet_g, 1, MPI_INT, MPI_SUM, &
                  0, COMM_model, MPIerr)
             call MPI_Reduce (nwet3d, nwet3d_g, 1, MPI_INT, MPI_SUM, &
                  0, COMM_model, MPIerr)

             if (mype_model==0) then
                write(*,'(a,5x,a,3x,i11)') &
                     'NEMO-PDAF', 'Number of global wet surface points', nwet_g
                write(*,'(a,5x,a,8x,i11)') &
                     'NEMO-PDAF', 'Number of global 3D wet points', nwet3d_g
                write(*,'(a,5x,a,8x,i11)') &
                     'NEMO-PDAF', 'global 2D wet points * nlayers', nwet_g*nk_p
             end if
          end if
       elseif (screen>2) then
          write(*,'(a,2x,a,1x,i4,2x,a,3x,i11)') &
               'NEMO-PDAF', 'PE', mype_model, 'Number of wet surface points', nwet
          write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
               'NEMO-PDAF', 'PE', mype_model, 'Number of 3D wet points', nwet3d
          write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
               'NEMO-PDAF', 'PE', mype_model, '2D wet points * nlayers', nwet*nk_p
       end if
    end if


! ******************************************************************
! *** Specify domain limits to limit observations to sub-domains ***
! ******************************************************************

    if (type_limcoords==0) then
       lim_coords(1,1) = glamt(i0 + 1, j0 + 1) * deg2rad
       lim_coords(1,2) = glamt(i0 + ni_p, j0 + 1) * deg2rad
       lim_coords(2,1) = gphit(i0 + ni_p, j0 + nj_p) * deg2rad
       lim_coords(2,2) = gphit(i0 + 1, j0 + 1) * deg2rad
    elseif (type_limcoords==1) then
       lim_coords(1,1) = minval(lon1_p) * deg2rad
       lim_coords(1,2) = maxval(lon1_p) * deg2rad
       lim_coords(2,1) = maxval(lat1_p) * deg2rad
       lim_coords(2,2) = minval(lat1_p) * deg2rad
    else
       lim_coords(1,1) = minval(glamt(:, :)) * deg2rad
       lim_coords(1,2) = maxval(glamt(:, :)) * deg2rad
       lim_coords(2,1) = maxval(gphit(:, :)) * deg2rad
       lim_coords(2,2) = minval(gphit(:, :)) * deg2rad
    end if

    call PDAFomi_set_domain_limits(lim_coords)

  end subroutine set_nemo_grid


#ifdef PDAF_OFFLINE
!> Return date as real from step value
!!
!! This routine mimicks what NEMO's routine calc_date would
!! return. However, here we assume that 'step' is just an integer
!! containing the date
!! 
  subroutine calc_date(step, rdate)

    implicit none

    integer, intent(in)    :: step
    real(pwp), intent(out) :: rdate

    rdate = REAL(step+1, pwp)
 
  end subroutine calc_date
#endif

end module nemo_pdaf
