!> Ensemble Initialization
!!
!! This routine calls the routines for initializing the ensemble.
!! 
!! The routine is called when the filter is initialized in
!! `PDAF_init`.
!! 
!! The routine is called by all filter processes and
!! initializes the ensemble for the process-local domain.
!! 
!!  **Calling Sequence**
!! 
!!  - Called from: `init_pdaf/PDAF_init` (PDAF module)
!!
subroutine init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  use mod_kind_pdaf
  use parallel_pdaf, &
       only: mype_filter
  use assimilation_pdaf, &
       only: dim_state, type_ens_init, type_central_state, ensscale, &
       coupling_nemo, screen
  use io_pdaf, &
       only: path_inistate, path_ens, file_ens, file_covar, &
       read_state_mv, read_ens_mv_loop, read_ens_states, read_ens_mv_filelist
  use statevector_pdaf, &
       only: n_fields, sfields
  use transforms_pdaf, &
       only: transform_field_mv
  use sample_ens_pdaf, &
       only: sample_ens_from_covar

  implicit none

! *** Arguments ***
  integer, intent(in) :: filtertype                     !< Type of filter to initialize
  integer, intent(in) :: dim_p                          !< PE-local state dimension
  integer, intent(in) :: dim_ens                        !< Size of ensemble
  real(pwp), intent(inout) :: state_p(dim_p)            !< PE-local model state
  !< It is not necessary to initialize the array 'state_p' for LETKF/LESTKF. 
  !< It is available here only for convenience and can be used freely.
  real(pwp), intent(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  real(pwp), intent(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  integer, intent(inout) :: flag                        !< PDAF status flag

! *** Local variables ***
  integer :: i, member              ! Counters
  integer :: id_field               ! Id of field in state vector
  real(pwp) :: ens_mean             ! Ensemble mean value
  real(pwp) :: inv_dim_ens          ! Inverse ensemble size
  integer :: verbose                ! Control verbosity


! ************************************
! *** Generate ensemble from files ***
! ************************************

  if (type_ens_init == 0) then

     ! Read ensemble states as model snapshots from a single file

     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble from file holding model snapshots'

     call read_ens_mv_loop(path_ens, dim_p, dim_ens, coupling_nemo, ens_p)

  elseif (type_ens_init == 1) then
     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble from single ensemble file'

     call read_ens_states(trim(file_ens), dim_p, dim_ens, ens_p)

  elseif (type_ens_init == 2) then
     
     ! Read ensemble states as model snapshots from separate files

     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble from list of output files'

      call read_ens_mv_filelist(1.0, .true., path_ens, dim_p, dim_ens, ens_p)

  elseif (type_ens_init == 3) then
     
     ! Initialize from covariance matrix

     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble by sampling from covariance matrix'

     state_p = 0.0
     call sample_ens_from_covar(trim(file_covar), dim_p, dim_ens, state_p, ens_p)

  elseif (type_ens_init == 4) then
     
     ! Ensemble restart using restarts read by NEMO

     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Ensemble restart'

     ! There is nothing to do in case of an ensemble restart

     ens_p = 0.0

     ! Deactivate replacing central state
     type_central_state = 0

  end if

  ! *** Options to replace the central state of the ensemble (i.e. the ensemble mean) ***

  if (type_central_state == 1) then

     ! Read ensemble central state vector state_p
     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Read central model state of ensemble from file'

     call read_state_mv(path_inistate, dim_p, 1, coupling_nemo, state_p)

  elseif (type_central_state == 2) then

     ! Obtain central state from model task 1 (set by model initialization)
     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Collect central ensemble state from model'

     call collect_state_init_pdaf(dim_p, state_p)

  end if


  ! *** Replace the ensemble central state if specified and scale variance of model fields ***

  if (type_central_state==1 .or. type_central_state==2) then

     ! *** Replace ensemble mean and inflate ensemble perturbations ***

     inv_dim_ens = 1.0_pwp/real(dim_ens, kind=pwp)

     if (mype_filter==0) write (*,'(a, 1x,a)') 'NEMO-PDAF', 'Set central state of ensemble'


     ! Scale ensemble perturbations - either all using 'ensscale' or field-wise
     do id_field = 1, n_fields

        ! Check general setting of ensscale
        if (ensscale /= 1.0 .and. sfields(id_field)%ensscale == 1.0) &
             sfields(id_field)%ensscale = ensscale

        if (mype_filter==0) &
             write (*,'(a, 1x,a,1x,a,1x,a,f12.4)') 'NEMO-PDAF', 'Scale field variance', &
             sfields(id_field)%variable, 'by', sfields(id_field)%ensscale

!$OMP PARALLEL DO private(member, ens_mean)
        do i = 1 + sfields(id_field)%off, sfields(id_field)%off + sfields(id_field)%dim
           ens_mean = 0.0_pwp
           do member = 1, dim_ens
              ens_mean = ens_mean + inv_dim_ens*ens_p(i, member)
           end do
        
           do member = 1, dim_ens
              ens_p(i, member) = sfields(id_field)%ensscale * (ens_p(i, member) - ens_mean) + state_p(i)
           end do
        end do
!$OMP END PARALLEL DO
     end do
     
  end if

  ! *** Transform fields

  do member = 1 , dim_ens

     if (mype_filter==0 .and. member==1) then
        verbose = 1
     else
        verbose = 0
     end if

     call transform_field_mv(1, ens_p(:,member), 11, verbose)
  end do

end subroutine init_ens_pdaf
