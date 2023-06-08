!> Initialize model dimensions
!!
!! Routine to perform initialization of the model information for
!! PDAF. Here, the global size of the model domain, it decomposition
!! and the coordinate arrays need to be initialized.
!!
subroutine init_grid_dims()

  use mod_kind_pdaf
  use mpi
  use netcdf
  use parallel_pdaf, &
       only: mype=>mype_world, npes=>npes_world, MPIerr, &
       abort_parallel
  use assimilation_pdaf, &
       only: dim_ens
  use nemo_pdaf, &
       only: path_dims, file_dims, tmask, jptra, &
       jpiglo, jpjglo, jpk, glamt, gphit, gdept_1d, &
       tmp_4d, dim_2d, dim_2d_p, lat1, lon1, &
       use_wet_state, nimpp, njmpp, istart, jstart, &
       nldi, nldj, nlei, nlej, ni_p, nj_p, nk_p, jptra
#if defined key_top
  use nemo_pdaf, &
       only: sn_tracer
#endif
  use io_pdaf, &
       only: verbose_io, check, add_slash
  use mod_memcount_pdaf, &
       only: memcount

  implicit none

! *** Arguments ***
  ! none

! *** local variables *** 
  integer :: screen=1                           ! Verbosity flag
  integer :: i, j, k                            ! Counters
  integer :: cnt, cnt_all, cnt_layers           ! Counters
  integer :: ncid                               ! nc file ID
  integer :: lon_dimid, lat_dimid, lvl_dimid    ! Dimension IDs
  integer :: id_gphit, id_glamt, id_navlev, id_temp  ! variable IDs
  real(8) :: missing_value                           ! missing value
  logical :: have_pdafnml                            ! Flag whether namelist file is present
  character(len=20) :: varname                  ! Name of variable ot read to determine mask
  character(len=1)  :: decomp='y'               ! Direction of decomposition (x, y)
  logical :: read_decomp=.false.                ! Whether to read domain decomposition from file
  character(len=50) :: file_decomp='decomp.txt' ! Name of decomposition file
  integer :: n_domains_lon, n_domains_lat
  integer, allocatable :: dims_lat(:), dims_lon(:)
  integer, allocatable :: nldi_all(:), nldj_all(:) ! first inner index in i/j direction for all PEs
  integer, allocatable :: nlei_all(:), nlej_all(:) ! last inner index in i/j direction for all PEs
  integer :: nx, ny, nz                         ! Size of 3D grid

  namelist /pdaf_offline/ screen, path_dims, file_dims, varname, use_wet_state, &
       read_decomp, file_decomp, verbose_io, jptra, dim_ens


! **********************
! *** Initialization ***
! **********************

  varname = 'thetao'
  path_dims = './'
  file_dims = 'nemo_output.nc'

! *** Read namelist file for PDAF ***

  open (500,file='pdaf.nml')
  read (500,NML=pdaf_offline)
  close (500)

  call add_slash(path_dims)

#if defined key_top
  if (jptra>0) allocate(sn_tracer(jptra))
#endif

! *************************************
! *** Read dimensions of model grid ***
! *************************************

  if (mype==0) then
     write (*,'(a,2x,a)') 'NEMO-PDAF', '*** Reading NEMO dimensions ***'
     write (*,'(a,2x,a,a)') 'NEMO-PDAF', 'File: ',trim(path_dims)//trim(file_dims)
     write (*,'(a,2x,a,a)') 'NEMO-PDAF', 'Variable used for mask: ', trim(varname)
  end if

  ! Open the file
  call check( nf90_open(trim(path_dims)//trim(file_dims), nf90_nowrite, ncid) )
    
  ! Get the dimensions
  call check( nf90_inq_dimid(ncid, 'x', lon_dimid) )  
  call check(nf90_inquire_dimension(ncid,lon_dimid,len=jpiglo))
    
  call check( nf90_inq_dimid(ncid, 'y', lat_dimid) )  
  call check(nf90_inquire_dimension(ncid,lat_dimid,len=jpjglo))
  
  call check( nf90_inq_dimid(ncid,'deptht', lvl_dimid) )  
  call check(nf90_inquire_dimension(ncid,lvl_dimid,len=jpk))

  nx = jpiglo
  ny = jpjglo
  nz = jpk


! ***********************************
! *** Define domain-decomposition ***
! ***********************************

  allocate(nldi_all(0:npes-1))
  allocate(nldj_all(0:npes-1))
  allocate(nlei_all(0:npes-1))
  allocate(nlej_all(0:npes-1))


  if (npes==1 .or. .not.read_decomp) then
     ! Use pre-defined decomposition

     decomp='y'
     if (npes==1) then
        ! No decomposition
        nldi_all(0) = 1
        nldj_all(0) = 1
        nlei_all(0) = nx
        nlej_all(0) = ny
     elseif (npes==2) then
        if (decomp=='y') then
           ! Decomposition in y
           nldi_all(:) = 1
           nlei_all(:) = nx
           nldj_all(0) = 1
           nlej_all(0) = ny/2
           nldj_all(1) = ny/2 + 1
           nlej_all(1) = ny
        else
           ! Decomposition in x
           nldj_all(:) = 1
           nlej_all(:) = ny
           nldi_all(0) = 1
           nlei_all(0) = nx/2
           nldi_all(1) = nx/2 + 1
           nlei_all(1) = nx
        end if
     else
        write (*,*) 'No decomposition defined in the code'
        call abort_parallel()
     end if
  else

     ! Read decomposition file

     allocate(dims_lon(0:npes-1))
     allocate(dims_lat(0:npes-1))

     if (mype==0) write(*,'(/a,2x,a,a)') &
          'NEMO-PDAF','*** Read domain decomposition: ', trim(file_decomp)

     open(11,FILE=trim(file_decomp))
     read(11,*) n_domains_lon, n_domains_lat
     if (mype==0) write (*,'(a,3x,a,2i5)') 'NEMO-PDAF', 'Number of tasks lon/lat: ', n_domains_lon, n_domains_lat

     do k = 0, n_domains_lon*n_domains_lat-1
        read(11, *)  nldj_all(k), nlej_all(k), nldi_all(k), nlei_all(k), dims_lon(k), dims_lat(k)
     end do
     close(11)

     if (n_domains_lon*n_domains_lat /= npes) then
        write (*,*) 'ERROR: number of domains in file inconsistent with number of processes'
        call abort_parallel()
     end if

     deallocate(dims_lon, dims_lat)

  end if

  ! Set sizes for this MPI task
  nldj = nldj_all(mype)
  nldi = nldi_all(mype)
  nlei = nlei_all(mype)
  nlej = nlej_all(mype)

  ! Set local domain size
  ni_p = nlei - nldi + 1
  nj_p = nlej - nldj + 1
  nk_p = jpk

  ! Start indices for sub-domain without halo
  istart = nldi
  jstart = nldj

  ! NEMO MPI offsets
  nimpp = 1
  njmpp = 1

  ! Dimension of 2D surface box
  dim_2d = jpiglo*jpjglo
  dim_2d_p = ni_p*nj_p

  ! Screen output
  if (mype==0 .and. npes>1 .and. screen>0) then
     write (*,'(/a,3x,a)') 'NEMO-PDAF','Grid decomposition:' 
     write (*,'(a, 8x,a,2x,a,a,2x,a,a,1x,a,6(1x,a))') &
          'NEMO-PDAF','rank ', 'istart', '  iend', 'jstart', '  jend', '  idim', '  jdim'
     do i = 0, npes-1
        write (*,'(a,2x, a,i6,1x,2i7,2i7,2i7)') 'NEMO-PDAF', 'RANK',i, nldj_all(i), nlej_all(i), &
             nldi_all(i), nlei_all(i), nlei_all(i)-nldi_all(i)+1, nlej_all(i)-nldj_all(i)+1
     end do
  end if


! *************************************************
! *** Read coordinates and field to define mask ***
! *************************************************

  call check( nf90_inq_varid(ncid, 'nav_lat', id_gphit) )
  call check( nf90_inq_varid(ncid, 'nav_lon', id_glamt) )
  call check( nf90_inq_varid(ncid, 'deptht', id_navlev) )
  call check( nf90_inq_varid(ncid, trim(varname), id_temp) )

  call check( nf90_get_att(ncid, id_temp, 'missing_value', missing_value) )

  allocate(glamt(nx, ny), gphit(nx, ny))
  allocate(gdept_1d(nz))
  allocate(tmp_4d(ni_p, nj_p, nz, 1), tmask(ni_p, nj_p, nz))
  call memcount(1, 'r', 2*nx*ny + nz + ni_p*nj_p*nz + ni_p*nj_p*nz)

  call check( nf90_get_var(ncid, id_glamt, glamt(:,:), (/1, 1/), (/nx, ny/) ) )
  call check( nf90_get_var(ncid, id_gphit, gphit(:,:), (/1, 1/), (/nx, ny/) ) )
  call check( nf90_get_var(ncid, id_navlev, gdept_1d(:), (/1/), (/nz/) ) )

  call check( nf90_get_var(ncid, id_temp, tmp_4d(:,:,:,:), (/nldi, nldj, 1, 1/), &
       (/ni_p, nj_p, nz, 1/) ) )
  
  ! Close the file. 
  call check( nf90_close(ncid) )


! *** Initialize mask  ***

  do k = 1, nz
     do j = 1, nj_p
        do i = 1, ni_p 
           if (abs(tmp_4d(i,j,k,1)- missing_value) > 0.1) then
              tmask(i,j,k) = 1.0_pwp
           else
              tmask(i,j,k) = 0.0_pwp
           endif
        enddo
     enddo
  end do


! *** Screen output ***

  if (mype==0 .and. screen>0) then
     write (*,'(/a,5x,a)') 'NEMO-PDAF', '*** NEMO: grid dimensions ***' 
     write(*,'(a,3x,2(6x,a),9x,a)') 'NEMO-PDAF', 'jpiglo','jpjglo','jpk' 
     write(*,'(a,3x,3i12)') 'NEMO-PDAF', jpiglo, jpjglo, jpk
     write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Dimension of global 3D grid box', jpiglo*jpjglo*jpk
     write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Number of global surface points', dim_2d
  end if

  if (npes>1 .and. screen>1) then
     write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
          'NEMO-PDAF', 'PE', mype, 'Dimension of local 3D grid box', ni_p*nj_p*nz
     write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
          'NEMO-PDAF', 'PE', mype, 'Number of local surface points', dim_2d_p
  end if

end subroutine init_grid_dims
