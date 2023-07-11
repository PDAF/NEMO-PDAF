module initialize_offline

  implicit none

  character(len=20) :: varname                  ! Name of variable to read to determine mask
  logical :: read_decomp=.false.                ! Whether to read domain decomposition from file
  character(len=50) :: file_decomp='decomp.txt' ! Name of decomposition file
  integer :: screen=1                           ! Verbosity flag

contains

!> Initialize parameters for PDAF offline implementation
!!
!! This routine reads in the pdaf offline namelist to 
!! initialize parameters for the offline implementation.
!! The routine afterwards calls the routine that initializes
!! the model grid information.
!!
  subroutine initialize

  use mod_kind_pdaf
  use parallel_pdaf, &
       only: mype_ens
  use assimilation_pdaf, &
       only: step_null
  use nemo_pdaf, &
       only: path_dims, file_dims, jptra, ndastp, use_wet_state, &
       type_limcoords
#if defined key_top
  use nemo_pdaf, &
       only: sn_tracer
#endif
  use io_pdaf, &
       only: verbose_io, check, add_slash, coupling_nemo, &
       incrTime, startIncrTime, endIncrTime

    implicit none

! *** Local variables ***    

    integer   :: day, year, month
    real(pwp) :: rdate

    namelist /pdaf_offline/ screen, path_dims, file_dims, varname, use_wet_state, &
         read_decomp, file_decomp, verbose_io, jptra, ndastp, &
         incrTime, startIncrTime, endIncrTime, coupling_nemo


! *** Initialize
    varname = 'thetao'
    path_dims = './'
    file_dims = 'nemo_grid.nc'

! *** Read namelist file for PDAF-offline ***

    open (500,file='pdaf_offline.nml')
    read (500,NML=pdaf_offline)
    close (500)

    call add_slash(path_dims)

    ! Set step_null
    rdate = real(ndastp)
    year = floor(rdate/10000.0_pwp)
    month = floor((rdate-real(year*10000))/100.0_pwp)
    day = floor(rdate-real(year*10000)-real(month*100))
    step_null = day

#if defined key_top
    if (jptra>0) allocate(sn_tracer(jptra))
#endif

    ! Print PDAF parameters to screen
    showconf: if (mype_ens == 0) then
       write (*, '(/a,1x,a)') 'NEMO-PDAF','-- Overview of PDAF-offline configuration --'
       write (*, '(a,3x,a)') 'NEMO-PDAF','[pdaf_offline]:'
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','screen       ', screen
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','use_wet_state', use_wet_state
       write (*, '(a,5x,a,6x,l)') 'NEMO-PDAF','read_decomp  ', read_decomp
       if (read_decomp) &
            write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','file_decomp  ', trim(file_decomp)
#if defined key_top
       write (*, '(a,5x,a,i10)') 'NEMO-PDAF','jptra        ', jptra
#endif
       write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','coupling_nemo', trim(coupling_nemo)
       write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','path_dims    ', trim(path_dims)
       write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','file_dims    ', trim(file_dims)
       write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','varname      ', trim(varname)
       write (*, '(a,5x,a,1x,i10)') 'NEMO-PDAF','ndastp       ', ndastp
       write (*, '(a,5x,a,f12.2)') 'NEMO-PDAF','incrTime       ', incrTime
       write (*, '(a,5x,a,f12.2)') 'NEMO-PDAF','startIncrTime  ', startincrTime
       write (*, '(a,5x,a,f12.2)') 'NEMO-PDAF','endIncrTime    ', endincrTime
       write (*, '(a,1x,a/)') 'NEMO-PDAF','-- End of PDAF-offline configuration overview --'
    end if showconf


! *** Initialize model grid information ***
    call init_grid_dims()

! *** Set initialization of lim_coords to using min/max of glamt/gphit
    type_limcoords = 2

  end subroutine initialize



!> Initialize model dimensions
!!
!! Routine to perform initialization of the model grid information
!! for PDAF. Here, the global size of the model domain, its
!! decomposition and the coordinate arrays are initialized.
!! 
!! This routine is particular for the offline implementation of
!! PDAF since the model grid information is not directly accessible
!! from the model.
!!
  subroutine init_grid_dims()

    use mod_kind_pdaf
    use mpi
    use netcdf
    use parallel_pdaf, &
         only: mype=>mype_world, npes=>npes_world, MPIerr, &
         abort_parallel
    use assimilation_pdaf, &
         only: step_null
    use nemo_pdaf, &
         only: path_dims, file_dims, tmask, tmp_4d, &
         jpiglo, jpjglo, jpk, glamt, gphit, gdept_1d, &
         dim_2d, dim_2d_p, nimpp, njmpp, istart, jstart, &
         nldi, nldj, nlei, nlej, ni_p, nj_p, nk_p, &
         lon1_p, lat1_p
    use io_pdaf, &
         only: check
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
    integer :: id_gphit, id_glamt, id_nlev, id_temp  ! variable IDs
    real(pwp) :: missing_value                    ! missing value
    integer :: n_domains_lon, n_domains_lat       ! Number of paralle tasks  in lon/lat
    integer, allocatable :: dims_lat(:), dims_lon(:)
    integer, allocatable :: nldi_all(:), nldj_all(:) ! first inner index in i/j direction for all PEs
    integer, allocatable :: nlei_all(:), nlej_all(:) ! last inner index in i/j direction for all PEs
    integer :: nx, ny, nz                         ! Size of 3D grid
    character(len=1)  :: decomp='y'               ! Direction of decomposition (x, y)
    real(pwp), allocatable :: glamt_g(:,:)        ! Global array for longitude coordinates
    real(pwp), allocatable :: gphit_g(:,:)        ! Global array for latitude coordinates
    real(pwp), allocatable :: lat1(:), lon1(:)    ! Coordinate vectors for regular lat/lon grids


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
       elseif (npes==4) then
          if (decomp=='b') then
             ! Decomposition in y
             nldi_all(0) = 1
             nlei_all(0) = nx/2
             nldi_all(1) = 1
             nlei_all(1) = nx/2
             nldi_all(2) = nx/2 + 1
             nlei_all(2) = nx
             nldi_all(3) = nx/2 + 1
             nlei_all(3) = nx
             nldj_all(0) = 1
             nlej_all(0) = ny/2
             nldj_all(1) = ny/2 + 1
             nlej_all(1) = ny
             nldj_all(2) = 1
             nlej_all(2) = ny/2
             nldj_all(3) = ny/2 + 1
             nlej_all(3) = ny
          elseif (decomp=='y') then
             ! Decomposition in y
             nldi_all(:) = 1
             nlei_all(:) = nx
             nldj_all(0) = 1
             nlej_all(0) = ny/4
             nldj_all(1) = ny/4 + 1
             nlej_all(1) = ny/2
             nldj_all(2) = ny/2 + 1
             nlej_all(2) = 3*ny/4
             nldj_all(3) = 3*ny/4 + 1
             nlej_all(3) = ny
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

    ! Set grid limits for this MPI task
    nldi = 1
    nldj = 1
    nlei = nlei_all(mype) - nldi_all(mype) + 1
    nlej = nlej_all(mype) - nldj_all(mype) + 1

    ! Set local domain size
    ni_p = nlei 
    nj_p = nlej 
    nk_p = jpk

    ! Start indices for sub-domain without halo
    istart = nldi_all(mype)
    jstart = nldj_all(mype)

    ! NEMO parallel grid offsets (with halo, but there is no halo in offline mode)
    nimpp = nldi_all(mype)
    njmpp = nldj_all(mype)

    ! Dimension of 2D surface box
    dim_2d = jpiglo*jpjglo
    dim_2d_p = ni_p*nj_p

    ! Screen output
    if (mype==0 .and. npes>1 .and. screen>0) then
       write (*,'(/a,3x,a)') 'NEMO-PDAF','Grid decomposition:' 
       write (*,'(a, 8x,a,2x,a,a,2x,a,a,1x,a,6(1x,a))') &
            'NEMO-PDAF','rank ', 'istart', '  iend', 'jstart', '  jend', '  idim', '  jdim'
       do i = 0, npes-1
          write (*,'(a,2x, a,i6,1x,2i7,2i7,2i7)') 'NEMO-PDAF', 'RANK',i, nldi_all(i), nlei_all(i), &
               nldj_all(i), nlej_all(i), nlei_all(i)-nldi_all(i)+1, nlej_all(i)-nldj_all(i)+1
       end do
    end if


! *************************************************
! *** Read coordinates and field to define mask ***
! *************************************************

    call check( nf90_inq_varid(ncid, 'nav_lat', id_gphit) )
    call check( nf90_inq_varid(ncid, 'nav_lon', id_glamt) )
    call check( nf90_inq_varid(ncid, 'deptht', id_nlev) )
    call check( nf90_inq_varid(ncid, trim(varname), id_temp) )

    call check( nf90_get_att(ncid, id_temp, 'missing_value', missing_value) )

    allocate(glamt_g(nx, ny), gphit_g(nx, ny))
    allocate(gdept_1d(nz))
    allocate(tmp_4d(ni_p, nj_p, nz, 1))
    call memcount(1, 'r', 2*nx*ny + nz + ni_p*nj_p*nz)

    call check( nf90_get_var(ncid, id_glamt, glamt_g(:,:), (/1, 1/), (/nx, ny/) ) )
    call check( nf90_get_var(ncid, id_gphit, gphit_g(:,:), (/1, 1/), (/nx, ny/) ) )
    call check( nf90_get_var(ncid, id_nlev, gdept_1d(:), (/1/), (/nz/) ) )

    call check( nf90_get_var(ncid, id_temp, tmp_4d(:,:,:,:), (/istart, jstart, 1, 1/), &
         (/ni_p, nj_p, nz, 1/) ) )
  
    ! Close the file. 
    call check( nf90_close(ncid) )


! *** Set coordinate vectors for rectangular grids ***

    allocate(lat1(ny), lon1(nx))
    allocate(lat1_p(nj_p), lon1_p(ni_p))
    call memcount(1, 'r', nx + ny + ni_p + nj_p)

    ! Global vectors
    lon1(:) = 0.0
    do j = 1, ny
       do i = 1, nx
          if (abs(glamt_g(i,j)) > 0.00001) then
             lon1(i) = glamt_g(i,j)
          endif
       enddo
    enddo

    lat1(:) = 0.0
    do j = 1, ny
       do i = 1, nx
          if (gphit_g(i,j) > 0.00001) then
             lat1(j) = gphit_g(i,j)
          endif
       enddo
    enddo

    ! Local vectors
    lon1_p = lon1(istart: istart + ni_p)
    lat1_p = lat1(jstart: jstart + nj_p)

    deallocate(lat1, lon1)


! *** Set local coordinate arrays ***

    allocate(glamt(ni_p, nj_p), gphit(ni_p, nj_p))
    call memcount(1, 'r', 2*ni_p*nj_p)

    glamt(1:ni_p, 1:nj_p) = glamt_g(istart: istart+ni_p-1, jstart: jstart+nj_p-1) 
    gphit(1:ni_p, 1:nj_p) = gphit_g(istart: istart+ni_p-1, jstart: jstart+nj_p-1) 


! *** Initialize mask ***

    allocate(tmask(ni_p, nj_p, nz))
    call memcount(1, 'r', ni_p*nj_p*nz)

    tmask = 0.0_pwp
    do k = 1, nz
       do j = 1, nj_p
          do i = 1, ni_p 
             if (abs(tmp_4d(i,j,k,1)- missing_value) > 0.1) then
                tmask(i,j,k) = 1.0_pwp
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

end module initialize_offline
