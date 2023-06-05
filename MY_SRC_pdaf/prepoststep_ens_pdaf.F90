!> Controlling Pre- and Post-Processing of the PDAF output
!!
!!  - For global filters (e.g. SEIK), the routine is called
!! before the analysis and after the ensemble transformation.
!!  - For local filters (e.g. LSEIK), the routine is called
!! before and after the loop over all local analysis
!! domains.
!! 
!! The routine provides full access to the state
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep
!! operations can be performed here. For example
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by
!! computing the estimated variances.
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!! 
!! If a user considers to perform adjustments to the
!! estimates (e.g. for balances), this routine is
!! the right place for it.
!! 
!! **Calling Sequence**
!! 
!!  - Called by: `PDAF_get_state` (as U_prepoststep) `PDAF_X_update` (as U_prepoststep)
!! 
subroutine prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  use mpi
  use mod_kind_pdaf
  use assimilation_pdaf, &
        only: step_null, ens_restart
  use parallel_pdaf, &
       only: mype=>mype_filter, comm_filter, MPIerr
  use statevector_pdaf, &
       only: n_fields, id, sfields
  use io_pdaf, &
        only: save_state, save_var, save_ens_sngl, file_PDAF_state, file_PDAF_variance, &
        write_field_mv, write_field_sngl, ids_write
  use nemo_pdaf, &
       only: ndastp, calc_date

  implicit none

! *** Arguments ***
  integer, intent(in) :: step           !< Current time step (negative for call after forecast)
  integer, intent(in) :: dim_p          !< PE-local state dimension
  integer, intent(in) :: dim_ens        !< Size of state ensemble
  integer, intent(in) :: dim_ens_p      !< PE-local size of ensemble
  integer, intent(in) :: dim_obs_p      !< Dimension of observation vector
  real, intent(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  real, intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  real, intent(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  integer, intent(in) :: flag           !< PDAF status flag


! *** local variables ***
  integer :: i, j, member, iens        ! counters
  integer :: id_var                    ! Index of a field variable in state vector
  logical, save :: firsttime = .true.  ! Routine is called for first time?
  real :: invdim_ens                   ! Inverse ensemble size
  real :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  real, allocatable :: state_tmp(:)    ! temporary state vector; holds state variances or increment
  integer,save :: writestep_var=1      ! Time index for file output of variance
  integer,save :: writestep_state=1    ! Time index for file output of state
  integer :: nsteps                    ! Number of steps written into file
  character(len=3) :: forana           ! String indicating forecast or analysis
  character(len=3) :: ensstr           ! Ensemble ID as string
  character(len=12) :: ndastp_str      ! String for model date 
  character(len=200) :: titleState, titleVar   ! Strings for file titles
  integer, allocatable :: dimfield_p(:) ! Local field dimensions
  integer, allocatable :: dimfield(:)  ! Global field dimensions
  real, allocatable :: rmse_est_p(:)   ! PE-local estimated RMS errors (ensemble standard deviations)
  real, allocatable :: rmse_est(:)     ! Global estimated RMS errors (ensemble standard deviations)
  logical :: inirestart                ! Whether prepoststep is called first time with ensemble restart
  real :: dimfield_inv                 ! Inverse of a field dimension
real :: rdate  
  
! **********************
! *** INITIALIZATION ***
! **********************

  inirestart = .false.
  if (step-step_null==0) then
     if (mype==0) write (*,'(a, i7,3x,a)') 'NEMO-PDAF', step, 'Analyze initial state ensemble'
     forana = 'ini'
     if (ens_restart) inirestart = .true.
  else if (step>0) then
     if (mype==0) write (*,'(a, 5x,a)') 'NEMO-PDAF', 'Analyze assimilated state ensemble'
     forana = 'ana'
  else
     if (mype==0) write (*,'(a, 5x,a)') 'NEMO-PDAF', 'Analyze forecast state ensemble'
     forana = 'for'
  end if

  ! Allocate fields
  allocate(state_tmp(dim_p))

  ! Initialize numbers
  invdim_ens    = 1.0_8 / real(dim_ens,8)  
  invdim_ensm1  = 1.0_8 / real(dim_ens - 1,8)



! **************************************************************
! *** Perform prepoststep for ensemble filter.               ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! **************************************************************

  ! *** Compute mean state

  state_p = 0.0
  do member = 1, dim_ens
!$OMP PARALLEL DO  
     do i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     end do
  end do
!$OMP PARALLEL DO  
  do i = 1, dim_p
    state_p(i) = invdim_ens * state_p(i)
  enddo


  ! *** Compute sampled variances ***
  state_tmp(:) = 0.0
  
  do member = 1, dim_ens
!$OMP PARALLEL DO  
     do j = 1, dim_p
        state_tmp(j) = state_tmp(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     end do
  end do
!$OMP PARALLEL DO  
  do j = 1, dim_p
     state_tmp(j) = invdim_ensm1 * state_tmp(j)
  enddo
  

! ******************************************************************
! *** Compute ensemble standard deviation (estimated RMS errors) ***
! ******************************************************************

  if (.not. inirestart) then

     allocate(rmse_est_p(n_fields))
     allocate(dimfield_p(n_fields))
     allocate(rmse_est(n_fields))
     allocate(dimfield(n_fields))

     ! Total sum of variance per field per process
     rmse_est_p  = 0.0
     do j = 1, n_fields
        do i = 1+sfields(j)%off, sfields(j)%dim+sfields(j)%off
           rmse_est_p(j) = rmse_est_p(j) + state_tmp(i)
        enddo
     enddo

     ! Global sum of variance per field
     call MPI_Reduce (rmse_est_p, rmse_est, n_fields, MPI_DOUBLE_PRECISION, MPI_SUM, &
          0, COMM_filter, MPIerr)

     ! Get global field dimensions
     dimfield_p(:) = sfields(:)%dim
     call MPI_Reduce(dimfield_p, dimfield, n_fields, MPI_INTEGER, MPI_SUM, 0, COMM_filter, MPIerr)

     ! Total estimated ensemble standard deviation per field per process
     if (mype == 0) then
        do j = 1, n_fields
           dimfield_inv = 1.0 / real(dimfield(j), 8)

           rmse_est(j) = rmse_est(j) * dimfield_inv
        end do

        ! Get global STDDEV
        rmse_est = sqrt(rmse_est)
     end if

     ! Output ensemble standard deviations
     if (mype == 0) then
        write (*, '(a,6x,a)') 'NEMO-PDAF', 'Ensemble standard deviation (estimated RMS error)'
        do i = 1, n_fields
           write (*,'(a,4x,a8,4x,a10,2x,es12.4)') &
                'NEMO-PDAF', 'RMSE-'//forana, trim(sfields(i)%variable), rmse_est(i)
        end do
     end if

  end if


! *******************
! *** File output ***
! *******************

  ! Set time string
  WRITE(ndastp_str,'(I8.8)') ndastp
  call calc_date(ABS(step), rdate)
!if (mype==0) write (*,*) 'RDATE', step, rdate
!  write(ndastp_str,'(F12.3)') rdate

  ! *** Write variance into nc file ***

  writevar: if (trim(save_var)=='ana' .or. trim(save_var)=='fcst' .or. trim(save_var)=='both') then

     if (writestep_var==1) then

        if (mype == 0 .and. .not. inirestart) then
           if (forana=='for' .and. (trim(save_var)=='fcst' .or. trim(save_var)=='both')) then
              write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance before analysis step'
           elseif (forana=='ana' .and. (trim(save_var)=='ana' .or. trim(save_var)=='both')) then
              if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance after analysis'
           else if (forana=='ini' .and. (.not.ens_restart)) then
              if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write initial variance file'
           end if
        end if

        if (save_var=='both' .and. forana/='ini') then
           nsteps = 2
        else
           nsteps = 1
        end if

        titleVar='Ensemble variance'

        if (forana/='ini') then
           if (trim(save_var)=='both' .or. (trim(save_var)=='fcst' .and. forana=='for') &
                .or. (trim(save_var)=='ana' .and. forana=='ana')) then
              call write_field_mv(state_tmp, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'.nc', &
                   titleVar, 1.0, nsteps, writestep_var, 0)
           end if
        else if (.not. (ens_restart)) then
           call write_field_mv(state_tmp, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'_ini.nc', &
                titleVar, 1.0, nsteps, writestep_var, 0)
        end if
        if (forana/='ini' .and. trim(save_var)=='both') writestep_var = 2

     elseif (writestep_var>1 .and. (trim(save_var)=='ana' .or. trim(save_var)=='both')) then
        
        if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance after analysis step'

        call write_field_mv(state_tmp, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'.nc', &
             titleVar, 1.0, 2, writestep_var, 0)

        writestep_var = 1
     end if
  end if writevar


  ! *** Write state into nc file ***
  writestate: if (trim(save_state)/='none' .and. (.not. (ens_restart .and. forana=='ini')) ) then

     ! Store state in state_tmp to avoid changing state_p
     state_tmp = state_p

     titleState = 'Ensemble mean state'

     if (trim(save_state)=='both' .and. forana/='ini') then
        nsteps = 2
     else
        nsteps = 1
     end if

     ! Write state file for viewing
     if (forana=='for' .and. (trim(save_state)=='fcst' .or. trim(save_state)=='both')) then
        if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble mean before analysis step'
        writestep_state = 1
     elseif (forana=='ana' .and. (trim(save_state)=='ana' .or. trim(save_state)=='both')) then
        if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble mean after analysis step'

        if (trim(save_state)=='both') then
           writestep_state = 2
        else
           writestep_state = 1
        end if
     elseif (forana=='ini') then
        if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble mean at initial time'
        writestep_state = 1
     else
        ! No file writing
        writestep_state = 0
     end if

     ! Write forecast and analysis into the same file
     if (forana/='ini') then
        if (writestep_state>0) then
           call write_field_mv(state_tmp, trim(file_PDAF_state)//'_'//trim(ndastp_str)//'.nc', &
                titleState, 1.0, nsteps, writestep_state, 1)
        end IF
     else
        call write_field_mv(state_tmp, trim(file_PDAF_state)//'_'//trim(ndastp_str)//'_ini.nc', &
             titleState, 1.0, nsteps, writestep_state, 1)
     end if
  endif writestate


  ! *** Write chlorophyll ensemble into nc files ***
  savens: if (save_ens_sngl .and. (.not. (ens_restart .and. forana=='ini')) ) then

     do i = 1, n_fields

        id_var = ids_write(i)

        if (id_var>0 .and. id_var<n_fields) then

           do iens = 1, dim_ens

              write(ensstr,'(i3.3)') iens

              ! Store state in state_tmp to avoid changing state_p
              state_tmp = ens_p(:,iens)

              titleState = 'Ensemble state vector'

              ! Write state file for viewing
              if (iens==1) then
                 if (forana=='for') then
                    if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble state before analysis step'
                    writestep_state = 1
                 elseif (forana=='ana') then
                    if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble state after analysis step'
                    writestep_state = 2
                 else
                    if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble state at initial time'
                    writestep_state = 1
                 end if
              end if

              ! Write forecast and analysis into the same file
              if (forana/='ini') then
                 call write_field_sngl(state_tmp, &
                      trim(file_PDAF_state)//'_'//trim(sfields(id_var)%variable)//'_'//trim(ndastp_str)//'_'//ensstr//'.nc', &
                      titleState, 1.0, 2, writestep_state, 1, id_var)
              else
                 call write_field_sngl(state_tmp, &
                      trim(file_PDAF_state)//'_'//trim(sfields(id_var)%variable)//'_'//trim(ndastp_str)//'_'//ensstr//'_ini.nc', &
                      titleState, 1.0, 1, writestep_state, 1, id_var)
              end if
           end do
        end if
     end do
  end if savens


! ********************
! *** finishing up ***
! ********************

  deallocate(state_tmp) 

  if (.not. inirestart) then
     deallocate(rmse_est_p, dimfield_p)
     if (mype==0) deallocate(rmse_est, dimfield)
  end if

  firsttime = .false.

end subroutine prepoststep_ens_pdaf
