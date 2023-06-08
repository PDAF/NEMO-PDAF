!> Module holding IO operations for NEMO-PDAF
!!
!! This code bases in wide parts on the implementation
!! by Wibke Duesterhoeft-Wriggers, BSH, Germany for the
!! CMEMS Baltic Monitoring and forecasting center
!!
module io_pdaf

  use mod_kind_pdaf
  use mpi

  ! Include dimension information for model grid
  use nemo_pdaf, &
       only: nlvls=>jpk, nlats=>jpjglo, nlons=>jpiglo, &
       depths=>gdept_1d, lons, lats, i0, j0, &
       tmp_4d, ni_p, nj_p, nk_p, istart, jstart, &
       nimpp, njmpp, nlei, nlej, stmp_4d, tmask

  ! Include information on state vector
  use statevector_pdaf, &
       only: id, sfields, n_fields, n_fields_covar

  ! Include parallelization information
  use parallel_pdaf, &
       only: mype=>mype_ens, npes=>npes_ens, comm_filter, abort_parallel

  ! Include transformation routines
  use transforms_pdaf, &
       only: field2state, field2state_missval, state2field, transform_field_mv

  ! Include coupling flag
  use assimilation_pdaf, &
       only: coupling_nemo

  implicit none
  save

  integer :: verbose_io=0   ! Set verbosity of IO routines (0,1,2,3)

  ! Control of IO
  character(len=4) :: save_var='none'        ! Write variance at 'fcst', 'ana', 'both', or 'none'
  character(len=4) :: save_state='both'      ! Write variance at 'fcst', 'ana', 'both', or 'none'
  logical :: save_ens_states=.false.         ! Write a single file of ensmeble state vectors
  logical :: save_ens_fields=.false.         ! Write set of files holding ensemble fields
  logical :: save_ens_sngl=.false.           ! write set of files holding ensemble of selected field
  logical :: save_incr                       ! Write increment to file
  logical :: do_deflate=.false.              ! Deflate variables in NC files (this seems to fail for parallel nc)
  character(len=3) :: sgldbl_io='sgl'        ! Write PDAF output in single (sgl) or double (dbl) precision

  character(len=100) :: file_PDAF_state='state'       ! File name for outputs of ensemble mean state
  character(len=100) :: file_PDAF_incr='incr'         ! File name for increment
  character(len=100) :: file_PDAF_variance='variance' ! File name for ensemble variance
  character(len=200) :: path_inistate      ! Path to NEMO files
  character(len=200) :: path_ens           ! Path of ensemble file
  character(len=200) :: file_ens           ! File name of ensemble file
  character(len=200) :: path_covar         ! Path of file holding covariance matrix
  character(len=200) :: file_covar         ! File name of file holding covariance matrix
  character(len=200) :: path_restart       ! Path of restart file
  character(len=80)  :: file_restart       ! file name of restart dile

  integer :: ids_write(25)

   ! Temporary - from offline code
  real(pwp) :: startEnsTime=1.0_pwp, endEnsTime=1.0_pwp, incrTime=1.0_pwp

  ! NEMO output file
  integer(4)        :: ntimec=1

  ! Missing value in netcdf file
  real(pwp) :: missing_value

contains
! ===================================================================================

!> Read fields from NEMO file into a state vector
!!
  subroutine read_state_mv(path, dim_p, itime, coupling, state)

    use netcdf

    implicit none

    character(len = *), intent(in)    :: path          !< Path to file
    integer(4),         intent(in)    :: dim_p         !< PE-local state dimension
    integer(4),         intent(in)    :: itime         !< Time to read in file
    character(len = *), intent(in)    :: coupling      !< Type of NEMO coupling
    real(pwp),          intent(inout) :: state(dim_p)  !< State vector

    ! Local variables
    integer(4) :: i               ! Counters
    integer(4) :: varid           ! Variable ID
    integer(4) :: ncid            ! NC file id
    character(len=50) :: filename ! Full file name


    if (verbose_io>0 .and. mype==0) &
         write(*,'(a,4x,a,i8)') 'NEMO-PDAF', '*** Read model output at time step: ', itime

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize state
    state = 0.0_pwp

    do i = 1, n_fields

       filename = trim(sfields(i)%file_state)
       if (verbose_io>1 .and. mype==0) then
          write(*,'(a,2x,a)') 'NEMO-PDAF', trim(path)//trim(filename)
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )
       if (coupling/='rest') then
          call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )
       else
          missing_value=0.0_pwp
       endif

       ! Read variable
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       if (sfields(i)%ndims == 3) then
          call check( nf90_get_var(ncid, varid, tmp_4d, &
               start=(/istart, jstart, 1, itime/), count=(/ni_p, nj_p, nlvls, 1/)) )
       else
          call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
               start=(/istart, jstart, itime/), count=(/ni_p, nj_p, 1/)) )
       end if

       call check( nf90_close(ncid) )

       ! Convert field to state vector
       call field2state_missval(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

    end do

    if (verbose_io>2) then
       do i = 1, n_fields
          write(*,'(a, 1x, a, a10, 1x, a,5x, 2f12.6)') &
               'NEMO-PDAF', 'Min and max for ',trim(sfields(i)%variable),' :     ',              &
               minval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
               maxval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
       enddo
    end if

  end subroutine read_state_mv


!===============================================================================

!> Read an ensemble of model fields into the ensemble array
!!
  subroutine read_ens_mv_loop(path, dim_p, dim_ens, coupling, ens)

    use netcdf

    implicit none

! *** Arguments ***
    character(len = *), intent(in)   :: path                 !< Path of file
    integer(4),         intent(in)   :: dim_p                !< State dimension
    integer(4),         intent(in)   :: dim_ens              !< Ensemble size
    character(len = *), intent(in)   :: coupling             !< Type of NEMO coupling
    real(pwp),          intent(inout):: ens(dim_p, dim_ens)  !< Ensemble array

! *** Local variables ***
    integer(4) :: i, member        ! Counters
    integer(4) :: ncid             ! NC file ID
    integer(4) :: varid            ! Variable ID
    character(len=50) :: filename  ! Full file name

    if (verbose_io>0 .and. mype==0) &
         write(*,'(a,4x,a)') 'NEMO-PDAF','*** Ensemble: Read model snapshots'

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize ensemble
    ens = 0.0_pwp

    do i = 1, n_fields

       filename = trim(sfields(i)%file)
       if (verbose_io>1 .and. mype==0) then
          write(*,'(a,2x,a)') 'NEMO-PDAF', trim(path)//trim(filename)
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       !  Read field
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       ! Read missing value
       if (coupling/='rest') then
          call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )
       else
          missing_value=0.0_pwp
       endif

       do member = 1, dim_ens

          if (verbose_io>0 .and. mype==0 .and. i==1) &
               write (*,'(a,4x,a,i6)') 'NEMO-PDAF','--- read member', member

          if (sfields(i)%ndims == 3) then
             call check( nf90_get_var(ncid, varid, tmp_4d, &
                  start=(/istart, jstart, 1, member/), count=(/ni_p, nj_p, nlvls, 1/)) )
          else
             call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
                  start=(/istart, jstart, member/), count=(/ni_p, nj_p, 1/)) )
          end if

          ! Convert field to state vector
          call field2state_missval(tmp_4d, ens(:,member), sfields(i)%off, sfields(i)%ndims, missing_value)

       enddo

       call check( nf90_close(ncid) )

    end do

    if (verbose_io>2) then
       do i = 1, n_fields
          write(*,'(a, 1x, a, a10, 1x, a,1x, 2es13.6)') &
               'NEMO-PDAF','Ensemble min and max for ',trim(sfields(i)%variable),' :     ', &
               minval(ens(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:)), &
               maxval(ens(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:))
       enddo
    end if

  end subroutine read_ens_mv_loop


!===============================================================================

!!> Read ensemble as state vectors from ensemble file
!!
  subroutine read_ens_states(ensfile_fullname, dim_state, dim_ens, ens)

    use netcdf

    implicit none

    character(len=*), intent(in)    :: ensfile_fullname        !< Name and path of ensemble file
    integer(4),       intent(in)    :: dim_state               !< PE-local state dimension
    integer(4),       intent(in)    :: dim_ens                 !< Ensemble size
    real(pwp),        intent(inout) :: ens(dim_state, dim_ens) !< Ensemble array

    ! Local variables
    integer(4) :: i                     ! Counter
    integer(4) :: ncid                  ! NC file id
    integer(4) :: dim_state_file        ! state dimension in file
    integer(4) :: dim_ens_file          ! Ensemble size in file
    character(len=400) :: varstr        ! String describing variables in state vector
    character(len=400) :: varstr_file   ! String describing variables in state vector
    integer(4) :: dimstate_dimid, dimens_dimid, ens_varid


! *** Generate string describing the state vector ***
    varstr = ''
    do i = 1, n_fields
       if (i==1) then
          varstr = trim(sfields(i)%variable)
       else
          varstr = trim(varstr)//' '//trim(sfields(i)%variable)
       endif
    end do

! *** Read file

    if (verbose_io>0 .and. mype==0) then
       write(*,'(1x,a,a)') "--- Read ensemble file: ", trim(ensfile_fullname)
    end if

    ! Open the file
    call check( nf90_open(ensfile_fullname, nf90_nowrite, ncid) )

    ! Read the string describing the state vector
    call check( nf90_get_att(ncid, NF90_GLOBAL, "state_fields", varstr_file) )

    ! Check consistency of state vector setup
    if (trim(varstr) == trim(varstr_file)) then

       ! Get the dimensions
       call check( nf90_inq_dimid(ncid, 'dim_state', dimstate_dimid) )
       call check(nf90_inquire_dimension(ncid,dimstate_dimid,len=dim_state_file))

       call check( nf90_inq_dimid(ncid, 'dim_ens', dimens_dimid) )
       call check(nf90_inquire_dimension(ncid,dimens_dimid,len=dim_ens_file))

       ! Check consistency of state dimension
       if (dim_state_file == dim_state) then

          ! Check consistency of ensemble size
          if (dim_ens_file >= dim_ens) then

             !  Read ensemble
             call check( nf90_inq_varid(ncid, 'ensemble', ens_varid) )

             call check( nf90_get_var(ncid, ens_varid, ens, start=(/1,1/),count=(/dim_state,dim_ens/)) )

             if (dim_ens_file> dim_ens) &
                  write (*,*) 'Notice: Ensemble in file is larger than dim_ens'
          else
             write (*,'(1x,a)') 'ERROR: Ensemble in file is too small'
             write (*,'(1x,a)')  'Stopping program!'
             call abort_parallel()
          end if

       else
          write (*,'(1x,a)') 'ERROR: inconsistent state dimension'
          write (*,'(1x,a)')  'Stopping program!'
          call abort_parallel()
       end if

    else
       write (*,'(1x,a)') 'ERROR: inconsistent variables in state'
       write (*,'(1x,a)')  'Stopping program!'
       call abort_parallel()
    end if

    call check( nf90_close(ncid) )

  end subroutine read_ens_states


!================================================================================

!> Initialize ensemble array from a list of NEMO output files
!!
  subroutine read_ens_mv_filelist(flate, zeromean, inpath, dim_p, dim_ens, ens)

  implicit none

! *** Arguments ***
  real(pwp),          intent(in)    :: flate        !< inflation
  character(len=*),   intent(in)    :: inpath       !< Path to input files
  integer(4),         intent(in)    :: dim_p        !< State dimension
  integer(4),         intent(in)    :: dim_ens      !< Ensemble size
  logical,            intent(in)    :: zeromean     !< Remove mean value
  real(pwp),          intent(inout) :: ens(:, :)    !< Ensemble array

! *** Local variables ***
  integer(4)        :: i, k, iens          ! Counters
  integer           :: ios                 ! Flag for file reading
  character(len=50) :: filenames(n_fields) ! Filename
  real(pwp)         :: ens_mean            ! Ensemble mean
  real(pwp)         :: invsteps            ! Inverse of ensemble size


! *** Read ensemble from files ***

  if (verbose_io>0 .and. mype==0) &
       write(*,'(/a,1x,a)') 'NEMO-PDAF', '*** Generating ensemble from output files ***'

  do i = 1, n_fields
    open (unit=10 + i, file=trim(sfields(i)%file), iostat=ios)
    if (ios /= 0) write(*,*) 'Could not open file ', trim(sfields(i)%file)
  enddo

  iens=0

  ensloop: do

     do i = 1, n_fields
        read (10 + i,*,iostat=ios) filenames(i)
        if (ios/=0) exit ensloop
     enddo

     do k =1, ntimec
        iens = iens + 1
        if (verbose_io>0 .and. mype==0) write (*,'(a,1x,a,i8)') 'NEMO-PDAF', '--- Read ensemble member', iens
        call read_ens_mv(inpath, filenames, dim_p, k, ens(:,iens))
     enddo

     if (iens==dim_ens) exit ensloop

  enddo ensloop

  do i = 1, n_fields
     close(10 + i)
  enddo

  ! Check ensemble size
  if (iens<dim_ens) then
     write (*,'(/1x,a)') 'ERROR: Available files less than ensemble size!'
     write (*,'(1x,a)')  'Stopping program!'
     call abort_parallel()
  end if


! *** Subtract ensemble mean and inflate ensemble perturbations ***

  if (zeromean) then

     if (verbose_io>0 .and. mype==0) &
          write(*,'(a,1x,a)') 'NEMO-PDAF', '--- Subtract mean of ensemble snapshots'


     invsteps = 1.0_pwp/real(dim_ens, kind=pwp)


!$OMP PARALLEL DO private(k, ens_mean)
     do k=1,dim_p
        ens_mean = 0.0_pwp
        do i=1,dim_ens
           ens_mean = ens_mean + invsteps*ens(k,i)
        end do

        do i=1,dim_ens
           ens(k,i) = flate*(ens(k,i)-ens_mean)
        end do
     end do
!$OMP END PARALLEL DO

  end if

end subroutine read_ens_mv_filelist


!===============================================================================

!> Read a model field into the state vector of an ensemble array
!!
  subroutine read_ens_mv(path, filenames, dim_p, itime, state)

    use netcdf

    implicit none

! *** Arguments ***
    character(len = *), intent(in)   :: path         !< Path of file
    character(len = *), intent(in)   :: filenames(:) !< filenames of all fields
    integer(4),         intent(in)   :: dim_p        !< State dimension
    integer(4),         intent(in)   :: itime        !< Time in file to read
    real(pwp),          intent(inout):: state(dim_p) !< State vector

! *** Local variables ***
    integer(4) :: i                ! Counters
    integer(4) :: ncid             ! NC file ID
    integer(4) :: varid            ! Variable ID
    character(len=50) :: filename  ! Full file name

    if (verbose_io>0 .and. mype==0) &
         write(*,'(a,4x,a,i8)') 'NEMO-PDAF','*** Ensemble: Read model output at time step: ', itime

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize state
    state = 0.0_pwp

    do i = 1, n_fields

       filename = filenames(i)

       if (verbose_io>1 .and. mype==0) then
          write(*,'(a,2x,a)') 'NEMO-PDAF', trim(path)//trim(filename)
          write (*,'(a,i5,1x,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       !  Read field
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       ! Read missing value
       call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )

       if (sfields(i)%ndims == 3) then
          call check( nf90_get_var(ncid, varid, tmp_4d, &
               start=(/istart, jstart, 1, itime/), count=(/ni_p, nj_p, nlvls, 1/)) )
       else
          call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
               start=(/istart, jstart, itime/), count=(/ni_p, nj_p, 1/)) )
       end if

       call check( nf90_close(ncid) )


       ! Convert field to state vector
       call field2state_missval(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

    end do

    if (verbose_io>2) then
       do i = 1, n_fields
          write(*,*) 'Min and max for ',trim(sfields(i)%variable),' :     ',              &
               minval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
               maxval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
       enddo
    end if

  end subroutine read_ens_mv

!================================================================================

!> Read initial forecast error covariance from a file
!!
  subroutine read_eof_cov(filename_cov, dim_state, dim_p, rank, state_p, eofV, svals, readmean)

    use nemo_pdaf, only: dim_2d_p, dim_3d_p
    use netcdf

! *** Arguments ***
    character(*), intent(in) :: filename_cov      !< filename
    integer, intent(in)      :: dim_state         !< global dimension of state vector
    integer, intent(in)      :: dim_p             !< local PE dimension of state vector
    integer, intent(in)      :: rank              !< rank of the covariance matrix
    real(pwp), intent(inout) :: eofV(dim_p, rank) !< singular vectors
    real(pwp), intent(inout) :: svals(rank)       !< singular values
    real(pwp), intent(inout) :: state_p(dim_p)    !< mean state vector
    logical, intent(in)      :: readmean          !< Whether to also read the mean state from covar file

! *** Local variables ***
    integer                :: ncid, dimid              ! file and dimension id
    integer                :: dim_file, rank_file      ! dimension size
    integer                :: varid                    ! variable id
    integer                :: maxvar                   ! Counter to distinguish reading EOFs and mean state
    integer                :: i_rank, i_field, i_var   ! counter
    integer                :: n_rank                   ! last dimension of variables
    integer                :: n_fields_read            ! Number of fields to read
    real(pwp)              :: missing_value=1.e+20_pwp ! missing value in variables
    character(len=5)       :: vartypes(2)              ! variable types: svd and mean

    ! *************************************************
    ! *** Initialize initial state and covar matrix ***
    ! *************************************************

    ! initialize the subroutine for reading
    vartypes(1) = '_svd'
    vartypes(2) = '_mean'

    if (.not.readmean) then
       ! Only read singular vectors
       maxvar = 1
    else
       ! Also read mean state vector from file
       maxvar = 2
    end if

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Open nc file
    if (verbose_io>0 .and. mype==0) &
         WRITE(*, '(a,1x,a,a)') 'NEMO-PDAF', '--- Reading covariance information from ', TRIM(filename_cov)

    call check( NF90_OPEN(TRIM(filename_cov), NF90_NOWRITE, ncid) )
!    call check( NF90_OPEN_PARopen_par(trim(filename_cov), NF90_NOWRITE, comm_filter, MPI_INFO_NULL, ncid) )

    ! Read rank stored in file
    call check( NF90_INQ_DIMID(ncid, 'rank', dimid) )
    call check( NF90_Inquire_dimension(ncid, dimid, len=rank_file) )

    ! Check consistency of dimensions
    checkdim: IF (rank_file < rank) THEN
      ! *** Rank stored in file is smaller than requested EOF rank ***
      WRITE(*, '(a)') 'Error: Rank stored in file is smaller than requested EOF rank'
      call check( NF90_CLOSE(ncid) )
      call abort_parallel()
    END IF checkdim

    ! read singular values
    call check( NF90_INQ_VARID(ncid, 'sigma', varid) )
    call check( NF90_GET_VAR(ncid, varid, svals, start=[1], count=[rank]) )

    ! Read singular vectors and mean
    ! position of the domain in the state vector
    do i_var = 1, maxvar
      if (i_var == 1) then
        n_rank = rank
      else if (i_var == 2) then
        n_rank = 1
      endif

      ! Determine number of fields to read from covariance matrix
      if (n_fields_covar <= 0) THEN
         n_fields_read = n_fields
      else
         n_fields_read = n_fields_covar
         if (verbose_io>0 .and. mype==0) &
              WRITE(*, '(a,1x,a,i4)') 'NEMO-PDAF', '--- Number of fields to read', n_fields_read
      end if

      ! loop over all fields
      do i_field = 1, n_fields_read

        if (verbose_io>1 .and. mype==0) &
             WRITE(*, '(a,1x,a,1x,a)') 'NEMO-PDAF', '--- read field', &
             trim(sfields(i_field)%variable)//trim(vartypes(i_var))

        call check( NF90_INQ_VARID(ncid, &
                   trim(sfields(i_field)%variable)//trim(vartypes(i_var)), &
                   varid) &
                  )
        do i_rank = 1, n_rank
          ! read and convert state vector into field shape
          if (sfields(i_field)%ndims == 2) then
            call check( NF90_GET_VAR(ncid, varid, tmp_4d(:, :, 1, 1), &
                        start=(/istart, jstart, i_rank/), count=(/ni_p, nj_p, 1/) ) )
          else if (sfields(i_field)%ndims == 3) then
            call check( NF90_GET_VAR(ncid, varid, tmp_4d(:, :, :, 1), &
                        start=(/istart, jstart, 1, i_rank/), count=(/ni_p, nj_p, nlvls, 1/) ) )
          end if

          if (i_var == 1) then
            ! convert field into state vector without dry points
            call field2state_missval(tmp_4d, eofV(:, i_rank), sfields(i_field)%off, &
                             sfields(i_field)%ndims, missing_value)
          else if (i_var == 2) then
            call field2state_missval(tmp_4d, state_p, sfields(i_field)%off, &
                             sfields(i_field)%ndims, missing_value)
          endif
        enddo

      enddo

    enddo
    call check( NF90_CLOSE(ncid) )

  end subroutine read_eof_cov

!================================================================================

!> Write a covariance matrix into file
!!
  subroutine write_covar_mv(path, name, dim_p, rank, svdU, meanstate, svals)
    
    use netcdf
   
    implicit none

! *** Arguments ***
    character(len=*), intent(in) :: path         ! File path
    character(len=*), intent(in) :: name         ! File name
    integer(4),       intent(in) :: dim_p        ! State dimension
    integer(4),       intent(in) :: rank         ! Number of EOFs to write
    real(8),          intent(inout) :: svdU(:,:) ! Array of singular vectors
    real(8),          intent(in) :: meanstate(:) ! State vector
    real(8),          intent(in) :: svals(:)     ! Vector of singular values

! *** Local variables ***
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: cnt,i,j,k, member
    integer(4) :: dimid_rank, dimid_lvls, dimid_lat, dimid_lon, dimid_one, dimid_state
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_field, id_svals
    integer(4) :: dimids(4)
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    real(8)    :: fillval
    real(8)    :: timeField(1)
    character(len=200) :: filename


! *** Create file ***

!    filename = trim(path)//trim(name)

    if (verbose_io>0 .and. mype==0) then
       write (*,'(1x,a)') 'Write covariance matrix into file'
       write (*,'(1x,a,a)') 'Create file: ', trim(path)//trim(name)
       if (do_deflate) write (*,'(1x,a)') '--- Apply deflation'
    endif

    if (npes==1) then
       call check( NF90_CREATE(trim(path)//trim(name),NF90_NETCDF4,ncid))
    else
       call check( NF90_CREATE_PAR(trim(path)//trim(name), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
    end if
    call check( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', 'Covariance matrix'))
    
    ! define dimensions for NEMO-input file
    call check( NF90_DEF_DIM(ncid,'rank', rank, dimid_rank))
    call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
    call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
    call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
    call check( NF90_DEF_DIM(ncid, 'one', 1, dimid_one) )
    call check( NF90_DEF_DIM(ncid, 'dim_state', dim_p, dimid_state) )
   
    dimids_field(4)=dimid_rank
    dimids_field(3)=dimid_lvls
    dimids_field(2)=dimid_lat
    dimids_field(1)=dimid_lon
       
    ! define variables
    call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
    call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
    call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))

    if (do_deflate) then
       call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
    end if

    call check( NF90_DEF_VAR(ncid, 'sigma', NF90_FLOAT, dimids_field(4), id_svals))

    do i = 1, n_fields
       if (sfields(i)%ndims==3) then
          dimids_field(3)=dimid_lvls
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable)//'_svd', NF90_DOUBLE, dimids_field(1:4), id_field) )
       else
          dimids_field(3)=dimid_rank
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable)//'_svd', NF90_DOUBLE, dimids_field(1:3), id_field) )
       end if
       if (do_deflate) &
            call check( NF90_def_var_deflate(ncid, id_field, 0, 1, 1) )
       call check( nf90_put_att(ncid, id_field, "coordinates", "nav_lat nav_lon") )
       fillval = 1.0e20
       call check( nf90_put_att(ncid, id_field, "_FillValue", fillval) )
       call check( nf90_put_att(ncid, id_field, "missing_value", fillval) )
    end do
       
    ! End define mode
    call check( NF90_ENDDEF(ncid) )

    ! write coordinates
    startz(1)=1
    countz(1)=nlvls

    startC(1)=1
    startC(2)=1
    countC(1)=nlons
    countC(2)=nlats

    if (mype==0) then
       call check( nf90_put_var(ncid,id_lev,depths,startz,countz))
       call check( nf90_put_var(ncid,id_lon,lons,startC,countC))
       call check( nf90_put_var(ncid,id_lat,lats,startC,countC))
    end if

     ! *** Write singular values

    if (mype==0) then
       startz(1)=1
       countz(1)=rank
       call check( nf90_put_var(ncid, id_svals, svals(1:rank), startz, countz))
    end if

    ! *** Write singular vectors as fields

    do member = 1, rank

       ! Backwards transformation of state
       call transform_field_mv(2, svdU(:,member),0, verbose_io)

    enddo

    do i = 1, n_fields

       if (verbose_io>1 .and. mype==0) then
          write (*,'(a,i5,a,1x,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       do member = 1, rank

          ! Convert state vector to field
          tmp_4d = 1.0e20
          call state2field(svdU(:,member), tmp_4d, sfields(i)%off, sfields(i)%ndims, tmask)

          if (verbose_io>0 .and. mype==0 .and. i==1) &
               write (*,'(a,4x,a,i6)') 'NEMO-PDAF','--- write rank', member

          call check( nf90_inq_varid(ncid, trim(sfields(i)%variable)//'_svd', id_field) )

          ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
          startt(1) = istart
          countt(1) = ni_p
          startt(2) = jstart
          countt(2) = nj_p
          startt(3) = 1
          countt(3) = nlvls
          startt(4) = member
          countt(4) = 1 
    
          if (sfields(i)%ndims==3) then
             startt(3) = 1
             countt(3) = nlvls
             call check( nf90_put_var(ncid, id_field, tmp_4d, startt, countt))
          else
             startt(3) = member
             countt(3) = 1 

             call check( nf90_put_var(ncid, id_field, tmp_4d, startt(1:3), countt(1:3)))
          end if

       end do
    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_covar_mv


!================================================================================
!> Write an ensemble file holding the state vectors
!!
  subroutine write_state_ens(file, dim_state, dim_ens, ens)

    use netcdf

    implicit none

! *** Arguments ***
    character(len=*), intent(in):: file          !< File name
    integer(4),       intent(in):: dim_state     !< state dimension
    integer(4),       intent(in):: dim_ens       !< Ensemble size
    real(pwp),        intent(in):: ens(:,:)      !< Ensemble array

! *** Local variables ***
    integer(4) :: i          ! Counter
    integer(4) :: fileid     ! NC file id
    integer(4) :: dimids(2)  ! dimension ids
    integer(4) :: id_ens     ! variable id
    real(pwp)  :: fillval    ! fill value
    integer(4) :: startv(2),countv(2)  ! Arrays for writing
    character(len=400) :: varstr   ! String describing variables in state vector
    character(len=200) :: filestr  ! String for file name

! *** Generate string describing the state vector ***
    varstr = ''
    do i = 1, n_fields
       if (i==1) then
          varstr = trim(sfields(i)%variable)
       else
          varstr = trim(varstr)//' '//trim(sfields(i)%variable)
       endif
    end do

    if (npes==1) then
       filestr = trim(file)//'.nc'
    else
       filestr = trim(file)//'_'//trim(str(mype))//'.nc'
    end if

! *** Write ensemble of state vectors ***

    ! *** Open file and initialize dimensions and fields ***
    call check( NF90_CREATE(trim(filestr),NF90_NETCDF4,fileid) )
    call check( NF90_PUT_ATT(fileid,NF90_GLOBAL,'title', &
         'Ensemble matrix for NEMO') )
    call check( nf90_put_att(fileid, NF90_GLOBAL, "state_fields", trim(varstr)) )

    ! define dimensions
    call check( NF90_DEF_DIM(fileid,'dim_state',dim_state,dimids(1)) )
    call check( NF90_DEF_DIM(fileid,'dim_ens',dim_ens,dimids(2)) )

    ! define variables
    call check( NF90_DEF_VAR(fileid,'ensemble',NF90_DOUBLE,dimids(1:2),id_ens) )
    fillval = 0.0_pwp
    call check( nf90_put_att(fileid, id_ens, "_FillValue", fillval) )
    call check( nf90_put_att(fileid, id_ens, "missing_value", fillval) )
    call check( NF90_def_var_deflate(fileid,id_ens,0,1,1) )

    ! End define mode
    call check( NF90_ENDDEF(fileid) )

    do i=1,dim_ens
       startv(1) = 1
       startv(2) = i
       countv(1) = dim_state
       countv(2) = 1
       call check( nf90_put_var(fileid,id_ens,ens(1:dim_state,i), startv,countv) )
    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(fileid) )

  end subroutine write_state_ens


!================================================================================

!> Write ensemble as single files holding model fields
!!
  subroutine write_ens_files(path,file_ens,dim_ens,ens)

    use netcdf

    implicit none

! *** Arguments ***
    character(len=*), intent(in) :: path        !< Path of file
    character(len=*), intent(in) :: file_ens    !< Name stub of file
    integer(4),       intent(in) :: dim_ens     !< Ensemble size
    real(pwp),     intent(inout) :: ens(:, :)   !< Ensemble array

! *** Local variables ***
    integer(4)          :: i                ! Counter
    character(len=200)  :: file_ensemble   ! Full name of an ensemble size
    character(len=200)  :: titleEns         ! NC title of file
    real(pwp)           :: time             ! Time in file


! *** Write ensemble perturbation files

    time=10.0_pwp !TO DO: this is random, time has to be read in pdaf.nml and set here

    titleEns='Ensemble perturbation (ens-mean) for PDAF'

    do i=1,dim_ens

       file_ensemble=trim(path)//trim(file_ens)//'_'//trim(str(i))//'.nc'

       call write_field_mv(ens(:, i), file_ensemble, titleEns,time, 1, 1, 1)
    enddo

  end subroutine write_ens_files

!================================================================================

!> Write a state vector as model fields into a file
!!
  subroutine write_field_mv(state, filename, title, &
       attime, nsteps, step, transform)

    use netcdf

    implicit none

! *** Arguments ***
    real(pwp),        intent(inout) :: state(:)  ! State vector
    character(len=*), intent(in) :: filename     ! File name
    character(len=*), intent(in) :: title        ! File title
    real(pwp),        intent(in) :: attime       ! Time attribute
    integer(4),       intent(in) :: nsteps       ! Number of time steps stored in file
    integer(4),       intent(in) :: step         ! Time index to write at
    integer(4),       intent(in) :: transform    ! Whether to transform fields

! *** Local variables ***
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: i
    integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon, dimid_one
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_field
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    integer(4) :: nf_prec      ! Precision for netcdf output of model fields
    real(pwp)  :: fillval
    real(4)    :: sfillval
    real(pwp)  :: timeField(1)
    integer(4) :: verbose      ! Control verbosity

    timeField(1)=attime

    if (sgldbl_io=='dbl') then
       if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))
       nf_prec = NF90_DOUBLE
    else
       if (.not. allocated(stmp_4d)) allocate(stmp_4d(ni_p, nj_p, nk_p, 1))
       nf_prec = NF90_FLOAT
    end if

    if (sgldbl_io=='dbl') then
       fillval = 1.0e20_pwp
    else
       sfillval = 1.0e20
    end if

    if (step==1) then

! *** Create file ***

       if (verbose_io>0 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Create file: ', trim(filename)

       if (npes==1) then
          call check( NF90_CREATE(trim(filename),NF90_NETCDF4,ncid))
       else
          call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
       end if
       call check( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', trim(title)))

       ! define dimensions for NEMO-input file
       call check( NF90_DEF_DIM(ncid,'t', nsteps, dimid_time))
       call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
       call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
       call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
       call check( NF90_DEF_DIM(ncid, 'one', 1, dimid_one) )

       dimids_field(4)=dimid_time
       dimids_field(3)=dimid_lvls
       dimids_field(2)=dimid_lat
       dimids_field(1)=dimid_lon

       ! define variables
       call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
       call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
       call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
       call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
       if (do_deflate) then
          call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
       end if

       do i = 1, n_fields
          if (sfields(i)%ndims==3) then
             dimids_field(3)=dimid_lvls
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), nf_prec, dimids_field(1:4), id_field) )
          else
             dimids_field(3)=dimid_time
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), nf_prec, dimids_field(1:3), id_field) )
          end if
          if (do_deflate) &
               call check( NF90_def_var_deflate(ncid, id_field, 0, 1, 1) )
          call check( nf90_put_att(ncid, id_field, "coordinates", "nav_lat nav_lon") )
          if (sgldbl_io=='dbl') then
             call check( nf90_put_att(ncid, id_field, "_FillValue", fillval) )
             call check( nf90_put_att(ncid, id_field, "missing_value", fillval) )
          else
             call check( nf90_put_att(ncid, id_field, "_FillValue", sfillval) )
             call check( nf90_put_att(ncid, id_field, "missing_value", sfillval) )
          end if
       end do

       ! End define mode
       call check( NF90_ENDDEF(ncid) )

       ! write coordinates
       startz(1)=1
       countz(1)=nlvls

       startC(1) = nimpp
       countC(1) = nlei-i0
       startC(2) = njmpp
       countC(2) = nlej-j0

       call check( nf90_put_var(ncid, id_lon, lons, startC, countC))
       call check( nf90_put_var(ncid, id_lat, lats, startC, countC))

       if (mype==0) then
          call check( nf90_put_var(ncid,id_lev,depths,startz,countz))
       end if

    else
       if (verbose_io>0 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Open file: ', trim(filename)

       if (npes==1) then
          call check( nf90_open(trim(filename), NF90_WRITE, ncid) )
       else
          call check( nf90_open_par(trim(filename), NF90_WRITE, comm_filter, MPI_INFO_NULL, ncid) )
       end if

    end if


    ! *** Write fields

    call check( nf90_inq_varid(ncid, 'time', id_time) )
    call check( nf90_VAR_PAR_ACCESS(NCID, id_time, NF90_COLLECTIVE) )
!    call check( nf90_put_vara(ncid, id_time, timeField, start=(/step/), count=(/1/)))
    startt(1) = step
    countt(1) = 1
    call check( nf90_put_var(ncid, id_time, timeField, startt(1:1), countt(1:1)))

    ! Backwards transformation of state fields
    if (mype==0) then
       verbose = 1
    else
       verbose = 0
    end if
    if (transform==1) call transform_field_mv(2, state, 0, verbose)

    do i = 1, n_fields

       ! Convert state vector to field
       if (sgldbl_io=='dbl') then
          tmp_4d = fillval
          call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims, tmask)
       else
          stmp_4d = sfillval
          call state2field(state, stmp_4d, sfields(i)%off, sfields(i)%ndims, tmask)
       end if

       if (verbose_io>1 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', '--- write variable: ', trim(sfields(i)%variable)
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), id_field) )
!       call check( nf90_VAR_PAR_ACCESS(NCID, id_field, NF90_COLLECTIVE) )

       ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = step
       countt(4) = 1

       if (sfields(i)%ndims==3) then
          startt(3) = 1
          countt(3) = nlvls

          if (sgldbl_io=='dbl') then
             call check( nf90_put_var(ncid, id_field, tmp_4d, startt, countt))
          else
             call check( nf90_put_var(ncid, id_field, stmp_4d, startt, countt))
          end if
       else
          startt(3) = step
          countt(3) = 1

          if (sgldbl_io=='dbl') then
             call check( nf90_put_var(ncid, id_field, tmp_4d, startt(1:3), countt(1:3)))
          else
             call check( nf90_put_var(ncid, id_field, stmp_4d, startt(1:3), countt(1:3)))
          end if
       end if

    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_field_mv


!================================================================================

!> Write a field from the state vector as model field into a file
!!
  subroutine write_field_sngl(state, filename, title, &
       attime, nsteps, step, transform, ifield)

    use netcdf

    implicit none

! *** Arguments ***
    real(pwp),        intent(inout) :: state(:)  ! State vector
    character(len=*), intent(in) :: filename     ! File name
    character(len=*), intent(in) :: title        ! File title
    real(pwp),        intent(in) :: attime       ! Time attribute
    integer(4),       intent(in) :: nsteps       ! Number of time steps stored in file
    integer(4),       intent(in) :: step         ! Time index to write at
    integer(4),       intent(in) :: transform    ! Whether to transform fields
    integer(4),       intent(in) :: ifield       ! ID of field to write

! *** Local variables ***
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: i
    integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon, dimid_one
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_field
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    integer(4) :: nf_prec      ! Precision for netcdf output of model fields
    real(pwp)  :: fillval
    real(4)    :: sfillval
    real(pwp)  :: timeField(1)
    integer(4) :: verbose      ! Control verbosity

    timeField(1)=attime

    if (sgldbl_io=='dbl') then
       if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))
       nf_prec = NF90_DOUBLE
    else
       if (.not. allocated(stmp_4d)) allocate(stmp_4d(ni_p, nj_p, nk_p, 1))
       nf_prec = NF90_FLOAT
    end if

    if (step==1) then

! *** Create file ***

       if (verbose_io>0 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Create file: ', trim(filename)

       if (npes==1) then
          call check( NF90_CREATE(trim(filename),NF90_NETCDF4,ncid))
       else
          call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
       end if
       call check( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', trim(title)))

       ! define dimensions for NEMO-input file
       call check( NF90_DEF_DIM(ncid,'t', nsteps, dimid_time))
       call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
       call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
       call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
       call check( NF90_DEF_DIM(ncid, 'one', 1, dimid_one) )

       dimids_field(4)=dimid_time
       dimids_field(3)=dimid_lvls
       dimids_field(2)=dimid_lat
       dimids_field(1)=dimid_lon

       ! define variables
       call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
       call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
       call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
       call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
       if (do_deflate) then
          call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
       end if

       do i = ifield, ifield
          if (sfields(i)%ndims==3) then
             dimids_field(3)=dimid_lvls
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), nf_prec, dimids_field(1:4), id_field) )
          else
             dimids_field(3)=dimid_time
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), nf_prec, dimids_field(1:3), id_field) )
          end if
          if (do_deflate) &
               call check( NF90_def_var_deflate(ncid, id_field, 0, 1, 1) )
          call check( nf90_put_att(ncid, id_field, "coordinates", "nav_lat nav_lon") )
          if (sgldbl_io=='dbl') then
             fillval = 1.0e20_pwp
             call check( nf90_put_att(ncid, id_field, "_FillValue", fillval) )
             call check( nf90_put_att(ncid, id_field, "missing_value", fillval) )
          else
             sfillval = 1.0e20
             call check( nf90_put_att(ncid, id_field, "_FillValue", sfillval) )
             call check( nf90_put_att(ncid, id_field, "missing_value", sfillval) )
          end if
       end do

       ! End define mode
       call check( NF90_ENDDEF(ncid) )

       ! write coordinates
       startz(1)=1
       countz(1)=nlvls

       startC(1) = nimpp
       countC(1) = nlei-i0
       startC(2) = njmpp
       countC(2) = nlej-j0

       call check( nf90_put_var(ncid, id_lon, lons, startC, countC))
       call check( nf90_put_var(ncid, id_lat, lats, startC, countC))

       if (mype==0) then
          call check( nf90_put_var(ncid,id_lev,depths,startz,countz))
       end if

    else
       if (verbose_io>0 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Open file: ', trim(filename)

       if (npes==1) then
          call check( nf90_open(trim(filename), NF90_WRITE, ncid) )
       else
          call check( nf90_open_par(trim(filename), NF90_WRITE, comm_filter, MPI_INFO_NULL, ncid) )
       end if

    end if


    ! *** Write fields

    call check( nf90_inq_varid(ncid, 'time', id_time) )
    call check( nf90_VAR_PAR_ACCESS(NCID, id_time, NF90_COLLECTIVE) )
!    call check( nf90_put_vara(ncid, id_time, timeField, start=(/step/), count=(/1/)))
    startt(1) = step
    countt(1) = 1
    call check( nf90_put_var(ncid, id_time, timeField, startt(1:1), countt(1:1)))

    ! Backwards transformation of state fields
    if (mype==0) then
       verbose = 1
    else
       verbose = 0
    end if
    if (transform==1) call transform_field_mv(2, state, 0, verbose)

    do i = ifield, ifield

       ! Convert state vector to field
       if (sgldbl_io=='dbl') then
          tmp_4d = 1.0e20_pwp
          call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims, tmask)
       else
          stmp_4d = 1.0e20
          call state2field(state, stmp_4d, sfields(i)%off, sfields(i)%ndims, tmask)
       end if

       if (verbose_io>1 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', '--- write variable: ', trim(sfields(i)%variable)
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), id_field) )
!       call check( nf90_VAR_PAR_ACCESS(NCID, id_field, NF90_COLLECTIVE) )

       ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = step
       countt(4) = 1

       if (sfields(i)%ndims==3) then
          startt(3) = 1
          countt(3) = nlvls

          if (sgldbl_io=='dbl') then
             call check( nf90_put_var(ncid, id_field, tmp_4d, startt, countt))
          else
             call check( nf90_put_var(ncid, id_field, stmp_4d, startt, countt))
          end if
       else
          startt(3) = step
          countt(3) = 1

          if (sgldbl_io=='dbl') then
             call check( nf90_put_var(ncid, id_field, tmp_4d, startt(1:3), countt(1:3)))
          else
             call check( nf90_put_var(ncid, id_field, stmp_4d, startt(1:3), countt(1:3)))
          end if
       end if

    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_field_sngl

!================================================================================

!> Write increment file
!!
  subroutine write_increment(state, filename, id_field)

    use netcdf

    implicit none

! *** Arguments ***
    real(pwp),     intent(inout) :: state(:)    !< State vector
    character(len=*), intent(in) :: filename    !< File name
    integer(4),       intent(in) :: id_field    !< Id of field in stte vector

! *** Local variables ***
    integer(4) :: i
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon
    integer(4) :: id_dateb, id_datef
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_incr
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    real(pwp)  :: fillval
    real(pwp)  :: bgnTimeInterv(1), finTimeInterv(1), timeInIncr(1)

    ! NOTE: This routine currently only writes a single field from the state vector
    !       which is specified by id_field. this must be a 3D field

    !The time values are read in via namelist at the moment
    ! Fixme
    timeInIncr(1)=incrTime !time for direct initialisation in Nemo (time of restart file which is used for adding to increment file)
    bgnTimeInterv(1)=startEnsTime !Start date of interval on which increment is valid (later for time ramp initialisation of  increment)
    finTimeInterv(1)=endEnsTime !End date of interval on which increment is valid (later for time rap init of increment)!

    if (verbose_io>0 .and. mype==0) &
         write (*,'(8x,a)') '--- Write increment file'

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    if (npes==1) then
       call check( NF90_CREATE(trim(filename), NF90_NETCDF4, ncid))
    else
       call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
    end if

    call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'title', &
         'Increment matrix for NEMO-PDAF coupling'))

    ! define dimensions for NEMO-input file
    call check( NF90_DEF_DIM(ncid, 't', 1, dimid_time))
    call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
    call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
    call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )

    dimids_field(4)=dimid_time
    dimids_field(3)=dimid_lvls
    dimids_field(2)=dimid_lat
    dimids_field(1)=dimid_lon

    ! define variables
    call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
    call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
    call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
    if (do_deflate) then
       call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
    end if

    call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
    call check( NF90_DEF_VAR(ncid, 'z_inc_dateb', NF90_DOUBLE, id_dateb))
    call check( NF90_DEF_VAR(ncid, 'z_inc_datef', NF90_DOUBLE, id_datef))

    do i = 1, n_fields

       if (sfields(i)%ndims==3) then
          dimids_field(3)=dimid_lvls
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), NF90_DOUBLE, dimids_field(1:4), id_incr) )
       else
          dimids_field(3)=dimid_time
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), NF90_DOUBLE, dimids_field(1:3), id_incr) )
       end if

       if (do_deflate) &
            call check( NF90_def_var_deflate(ncid, id_incr, 0, 1, 1) )

       fillval = 0.0_pwp
       call check( nf90_put_att(ncid, id_incr, "long_name", trim(sfields(i)%name_incr)//trim('Increment')) )
       call check( nf90_put_att(ncid, id_incr, "units", trim(sfields(i)%unit)) )
       call check( nf90_put_att(ncid, id_incr, "coordinates", "nav_lat nav_lon") )
       call check( nf90_put_att(ncid, id_incr, "_FillValue", fillval) )
       call check( nf90_put_att(ncid, id_incr, "missing_value", fillval) )
    end do

    ! End define mode
    call check( NF90_ENDDEF(ncid) )


    ! *** write coordinates ***
    startz(1)=1
    countz(1)=nlvls

    startC(1)=1
    startC(2)=1
    countC(1)=nlons
    countC(2)=nlats

    if (mype==0) then
       call check( nf90_put_var(ncid, id_lev, depths, startz, countz))
       call check( nf90_put_var(ncid, id_lon, lons, startC, countC))
       call check( nf90_put_var(ncid, id_lat, lats, startC, countC))

       call check( nf90_put_var(ncid, id_time, timeInIncr, start=(/1/), count=(/1/)))
       call check( nf90_put_var(ncid, id_dateb, bgnTimeInterv, start=(/1/), count=(/1/)))
       call check( nf90_put_var(ncid, id_datef, finTimeInterv, start=(/1/), count=(/1/)))
    end if

    ! *** Write fields ***

    do i = 1, n_fields

       ! Convert state vector to field
       tmp_4d = 0.0_pwp
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (verbose_io>1) &
            write (*,'(5x,a,a)') '--- write variable: ', trim(sfields(i)%name_incr)
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_incr), id_incr) )

       ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = 1
       countt(4) = 1

       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_incr, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_put_var(ncid, id_incr, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_increment


!================================================================================

!> Overwrite the NEMO restart file
!!
  subroutine update_restart_mv(state, state_tmp)

    use netcdf

    implicit none

! *** Arguments ***
    real(pwp),     intent(inout) :: state(:)     !< State vector
    real(pwp),     intent(inout) :: state_tmp(:) !< tmp state vector (used for storage)

! *** Local variables ***
    integer :: i                      ! Counter
    integer :: ncid                   ! NC file ID
    integer :: lid, uid               ! index range in state vector
    integer :: id_var_n, id_var_b     ! NC variable IDs
    integer  :: startt(4), countt(4)  ! arrays for file writing
    character(len=30) :: rst_file     ! Name of restart file
    integer(4) :: verbose      ! Control verbosity


    ! Attention in run script copy restart file from time of DA to file 'restart_trc_in_befDA.nc'

    !Write oxy to TRNOXY of restart file (now, time t) ->
    !restart Nemo with nn_euler=0 (TRBOXY is oxy for t-Delta t)

    ! Store name of restart file
    rst_file = sfields(1)%rst_file

    if (verbose_io>0) &
         write (*,'(a,3x,a,1x,a)') 'NEMO-PDAF', '--- Overwrite restart file:',trim(path_restart)//trim(rst_file)

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Open file and retrieve field ids
    if (npes==1) then
       call check( nf90_open(trim(path_restart)//trim(rst_file),NF90_WRITE, ncid))
    else
       call check( nf90_open_par(trim(path_restart)//trim(rst_file),NF90_WRITE,comm_filter, MPI_INFO_NULL, ncid))
    end if

    ! field transformation
    if (mype==0) then
       verbose = 1
    else
       verbose = 0
    end if
    call transform_field_mv(2, state, 21, verbose)

    do i = 1, n_fields

       if (trim(sfields(i)%rst_file) /= trim(rst_file)) then
       ! Open other restart file and retrieve field ids
          if (verbose_io>0 .and. mype==0) &
               write (*,'(a, 3x,a,1x,a)') 'NEMO-PDAF', '--- Open restart file:',trim(path_restart)//trim(sfields(i)%rst_file)
          if (npes==1) then
             call check( nf90_open(trim(path_restart)//trim(sfields(i)%rst_file),NF90_WRITE, ncid))
          else
             call check( nf90_open_par(trim(path_restart)//trim(sfields(i)%rst_file), &
                  NF90_WRITE, comm_filter, MPI_INFO_NULL, ncid))
          end if

          ! Store name of restart file
          rst_file = sfields(i)%rst_file

       end if

       ! Retrieve field IDs
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_n), id_var_n))
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_b), id_var_b))


       ! backwards transformation state - only if not done by write_increment before

       ! Convert state vector to field
       tmp_4d = 0.0_pwp
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       ! *** write variable for current time ***
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = 1
       countt(4) = 1

       if (sfields(i)%ndims==3) then
          call check( nf90_put_var(ncid, id_var_n, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_put_var(ncid, id_var_n, tmp_4d, startt(1:3), countt(1:3)))
       end if

       ! *** For second (past) time use increment ***

       ! Read field, add increment, and write field
       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_get_var(ncid, id_var_b, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_get_var(ncid, id_var_b, tmp_4d, startt(1:3), countt(1:3)))
       end if

       call field2state(tmp_4d, state, sfields(i)%off, sfields(i)%ndims)

       lid = sfields(i)%off+1
       uid = sfields(i)%off+sfields(i)%dim
       state(lid : uid) = state(lid : uid) + state_tmp(lid : uid)

       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_var_b, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_put_var(ncid, id_var_b, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    call check( nf90_close(ncid))

  end subroutine update_restart_mv


!================================================================================

!> Check status of NC operation
!!
  subroutine check(status)

    use netcdf

! *** Aruments ***
    integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       call abort_parallel()
    end if

  end subroutine check


! ==============================================================================

!> Add a trailing slash to a path string
!!
!! This routine ensures that a string defining a path
!! has a trailing slash.
!!
  subroutine add_slash(path)

    implicit none

! *** Arguments ***
    character(len=100) :: path  !< String holding the path

! *** Local variables ***
    integer :: strlength

! *** Add trailing slash ***
    strlength = len_trim(path)

    if (path(strlength:strlength) /= '/') then
       path = trim(path) // '/'
    end if

  end subroutine add_slash


!===============================================================================

!> Convert an integer to a strong of length 4
!!
  character(len=4) function str(k)

    implicit none

    integer, intent(in) :: k   !< number

    write (str, '(i4.4)') k

  end function str


!===============================================================================

!> Check whether a file exists
!!
  function file_exists(filename) result(res)

    implicit none

    character(len=*),intent(in) :: filename   !< File name
    logical                     :: res        !< Status of file

    ! Check if the file exists
    inquire( file=trim(filename), exist=res )

  end function file_exists

end module io_pdaf
