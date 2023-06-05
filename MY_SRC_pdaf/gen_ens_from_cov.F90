module sample_ens_pdaf

!> Generate ensemble perturbations from covariance matrix
!!
!! This routine read the covariance matrix information from a
!! file and then calls PDAF_sampleens to generate ensemble perturbations.
!! 

contains

  subroutine sample_ens_from_covar(filename_cov, dim_p, dim_ens, state_p, ens_p)

    use mod_kind_pdaf
    use pdaf_interfaces_module, &
         only: PDAF_SampleEns
    use io_pdaf, &
         only: read_eof_cov
    use assimilation_pdaf, &
         only: screen, dim_state
    use parallel_pdaf, &
         only: mype=>mype_filter, abort_parallel

    implicit none 

! *** Arguments ***
    character(*), intent(in) :: filename_cov          !< covariance filename
    integer, intent(in)      :: dim_p                 !< dimension of local state vector
    integer, intent(in)      :: dim_ens               !< ensemble size
    real(pwp), intent(inout) :: state_p(dim_p)        !< state vector
    real(pwp), intent(inout) :: ens_p(dim_p, dim_ens) !< ensemble array

! *** Local variables ***
    integer :: rank                         ! Rank of variance matrix read from file
    integer :: status_pdaf                  ! PDAF status flag
    integer :: verbose_sampleens            ! Set verbosity of PDAF_sampleens
    real(pwp), allocatable :: eofV(:, :)    ! Array holding singular vectors
    real(pwp), allocatable :: svals(:)      ! Vector holding singular values
    logical :: readmean                     ! Whther to read the ensemble mean state from covariance file


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

    ! *** Rank of matrix is ensemble size minus one
    rank = dim_ens - 1

    ! allocate memory for temporary fields
    allocate(eofV(dim_p, rank))
    allocate(svals(rank))

    ! get eigenvalue and eigenvectors from file
    readmean = .false.
    eofV = 0.0
    state_p = 0.0
    call read_eof_cov(filename_cov, dim_state, dim_p, rank, state_p, eofV, svals, readmean)

    ! *** Generate full ensemble on filter-PE 0 ***
    verbose_sampleens = 0
    if (mype==0) then
       write (*, '(a, 1x, a)') 'NEMO-PDAF', '--- generate ensemble using PDAF_SampleEns'
       write (*, '(a, 3x, a)') &
            'NEMO-PDAF', '--- use 2nd order exact sampling'
       write (*, '(a, 3x, a, i5)') 'NEMO-PDAF', '--- Ensemble size:  ', dim_ens
       write (*, '(a, 3x, a, i5)') 'NEMO-PDAF', '--- number of EOFs: ', rank

       if (screen>0) verbose_sampleens = 1
    endif

    ! Use PDAF routine to generate ensemble from covariance matrix
    call PDAF_SampleEns(dim_p, dim_ens, eofV, svals, state_p, ens_p, verbose_sampleens, status_pdaf)

    if (status_pdaf /= 0) then
       write (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in sample ensemble of PDAF - stopping! (PE ', mype, ')'
       call abort_parallel()
    end if


! ****************
! *** clean up ***
! ****************

    deallocate(svals, eofV)

  end subroutine sample_ens_from_covar

end module sample_ens_pdaf
