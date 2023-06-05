!> Assimilation Parameters
!!
!! This module provides variables needed for the
!! assimilation.
!! 
!! See `initialize_pdaf` for where many of these
!! variables are initialised.
!! 
module assimilation_pdaf

  use mod_kind_pdaf

  implicit none
  save

  integer :: type_ens_init = 0    !< Type of ensemble initialization
       !< (0) read snapshots from a single model file
       !< (1) read states from single ensemble file,
       !< (2) read snapshots from separate model files
       !< (3) initialize ensemble from covariance matrix file
       !< (4) ensemble restart using fields from NEMO restart files
  integer :: type_central_state = 1    !< Type of central state of ensemble
       !< (0) mean of model snapshots
       !< (1) read from file
       !< (2) use collect_state
  real(pwp) :: ensscale=1.0            !< Scaling factor for initial ensemble
       !< Type of coupling between NEMO and PDAF
  character(len=4)   :: coupling_nemo = 'odir'   !< offline: 'rest', 'incr', online: 'oinc', 'odir'
  logical :: ens_restart = .false.     !< Whether to perform ensemble restart using NEMO's restart


! *** Model- and data specific variables ***

  integer :: dim_state     !< Global model state dimension
  integer :: dim_state_p   !< Model state dimension for PE-local domain

  ! Settings for time stepping - available as namelist read-in
  integer :: step_null = 0       !< initial time step of assimilation
  logical :: restart = .false.   !< Whether to restart the data assimilation from previous run

! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! Settings for time stepping - available as command line options
  logical :: model_error     !< Control application of model error
  real(pwp) :: model_err_amp !< Amplitude for model error

! Settings for observations - available as command line options
  integer :: delt_obs         !< time step interval between assimilation steps
  logical :: twin_experiment  !< Whether to run an twin experiment with synthetic observations

! General control of PDAF - available as command line options
  integer :: screen       !< Control verbosity of PDAF
                          !< * (0) no outputs
                          !< * (1) progress info
                          !< * (2) add timings
                          !< * (3) debugging output
  integer :: dim_ens      !< Size of ensemble
  integer :: filtertype   !< Select filter algorithm:
                          !<   * SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                          !<   PF (12), GENOBS (100), 3DVAR (200)
  integer :: subtype      !< Subtype of filter algorithm
                          !<   * SEEK: 
                          !<     (0) evolve normalized modes
                          !<     (1) evolve scaled modes with unit U
                          !<     (2) fixed basis (V); variable U matrix
                          !<     (3) fixed covar matrix (V,U kept static)
                          !<   * SEIK:
                          !<     (0) ensemble forecast; new formulation
                          !<     (1) ensemble forecast; old formulation
                          !<     (2) fixed error space basis
                          !<     (3) fixed state covariance matrix
                          !<     (4) SEIK with ensemble transformation
                          !<   * EnKF:
                          !<     (0) analysis for large observation dimension
                          !<     (1) analysis for small observation dimension
                          !<   * LSEIK:
                          !<     (0) ensemble forecast;
                          !<     (2) fixed error space basis
                          !<     (3) fixed state covariance matrix
                          !<     (4) LSEIK with ensemble transformation
                          !<   * ETKF:
                          !<     (0) ETKF using T-matrix like SEIK
                          !<     (1) ETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LETKF:
                          !<     (0) LETKF using T-matrix like SEIK
                          !<     (1) LETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * ESTKF:
                          !<     (0) standard ESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LESTKF:
                          !<     (0) standard LESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * NETF:
                          !<     (0) standard NETF 
                          !<   * LNETF:
                          !<     (0) standard LNETF 
                          !<   * PF:
                          !<     (0) standard PF 
                          !<   * 3D-Var:
                          !<     (0) parameterized 3D-Var
                          !<     (1) 3D Ensemble Var using LESTKF for ensemble update
                          !<     (4) 3D Ensemble Var using ESTKF for ensemble update
                          !<     (6) hybrid 3D-Var using LESTKF for ensemble update
                          !<     (7) hybrid 3D-Var using ESTKF for ensemble update
  integer :: incremental  !< Perform incremental updating in LSEIK
  integer :: dim_lag      !< Number of time instances for smoother

! Filter settings - available as command line options
!    ! General
  integer :: type_forget  !< Type of forgetting factor
  real(pwp) :: forget     !< Forgetting factor for filter analysis
  integer :: dim_bias     !< dimension of bias vector

!    ! ENKF
  integer :: rank_analysis_enkf  !< Rank to be considered for inversion of HPH

!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  integer :: type_trans    !< Type of ensemble transformation 
                           !< * SEIK/LSEIK: 
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ETKF/LETKF with subtype=4: 
                           !< (0) use deterministic symmetric transformation
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ESTKF/LESTKF:
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T
                           !< * NETF/LNETF:
                           !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           !< (1) use identity transformation

!    ! LSEIK/LETKF/LESTKF/LNETF
  integer :: locweight     !< Type of localizing weighting of observations
                    !<   * (0) constant weight of 1
                    !<   * (1) exponentially decreasing with SRANGE
                    !<   * (2) use 5th-order polynomial
                    !<   * (3) regulated localization of R with mean error variance
                    !<   * (4) regulated localization of R with single-point error variance
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  integer :: type_sqrt     !< Type of the transform matrix square-root 
                    !<   * (0) symmetric square root
                    !<   * (1) Cholesky decomposition
!    ! 3D-Var
  integer :: type_opt      !< Type of minimizer for 3DVar
                    !<   * (1) LBFGS (default)
                    !<   * (2) CG+
                    !<   * (3) plain CG
                    !<   * (12) CG+ parallelized
                    !<   * (13) plain CG parallelized
  integer :: dim_cvec = 0  !< Size of control vector (parameterized part; for subtypes 0,1)
  integer :: dim_cvec_ens = 0   !< Size of control vector (ensemble part; for subtypes 1,2)
  integer :: mcols_cvec_ens = 1 !< Multiplication factor for number of columns for ensemble control vector
  real(pwp) :: beta_3dvar = 0.5 !< Hybrid weight for hybrid 3D-Var
!    ! NETF/LNETF
  integer :: type_winf          !< Set weights inflation: (1) activate
  real(pwp)    :: limit_winf    !< Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  integer :: type_hyb      !< Type of hybrid weight: (2) adaptive from N_eff/N
  real    :: hyb_gamma     !< Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  real    :: hyb_kappa     !< Hybrid norm for using skewness and kurtosis
!    ! Particle filter
  integer :: pf_res_type   !< Resampling type for PF
                           !< (1) probabilistic resampling
                           !< (2) stochastic universal resampling
                           !< (3) residual resampling        
  integer :: pf_noise_type    !< Resampling type for PF
                           !< (0) no perturbations, (1) constant stddev, 
                           !< (2) amplitude of stddev relative of ensemble variance
  real(pwp) :: pf_noise_amp !< Noise amplitude (>=0.0, only used if pf_noise_type>0)

!    ! Other variables - _NOT_ available as command line options!
  real(pwp) :: time        !< model time

  integer, allocatable :: id_lstate_in_pstate(:) !< Indices of local state vector in global vector
  real(pwp) :: domain_coords(2)       !< Coordinates of local analysis domain

  integer :: assim_flag = 0   !< Flag whether assimilation step was just done

!$OMP THREADPRIVATE(domain_coords, id_lstate_in_pstate)

contains

  !> Performing the Assimilation Step
  !!
  !! This routine is called during the model integrations at each timestep.
  !! It calls PDAF to check whether the forecast phase is completed and if
  !! so, PDAF will perform the analysis step.
  !! 
  !! **Calling Sequence**
  !! 
  !!  - Called from: `step.F90'
  !! 
  !!  - Calls: `PDAFomi_assimilate_local`
  !! 
  subroutine assimilate_pdaf(kt)

    use pdaf_interfaces_module, &
         only: PDAFomi_assimilate_local, PDAF_get_localfilter
    use parallel_pdaf, &
         only: mype_ens, abort_parallel, COMM_ensemble, MPIerr
    use nemo_pdaf, &
         only: lwp, numout
    use asminc_pdaf, &
         only: update_asm_step_pdaf

    implicit none

! *** Arguments ***
    integer, intent(in) :: kt  ! time step

! *** Local variables ***
    integer :: status_pdaf         ! PDAF status flag
    integer :: localfilter         ! Flag for domain-localized filter (1=true)

    !! External subroutines 
    !!  (subroutine names are passed over to PDAF in the calls to 
    !!  PDAF_get_state and PDAF_assimilate_X. This allows the user 
    !!  to specify the actual name of a routine. However, the 
    !!  PDAF-internal name of a subroutine might be different from
    !!  the external name!)

    ! Interface between model and PDAF, and prepoststep
    external :: collect_state_pdaf, &  ! Collect a state vector from model fields
         distribute_state_pdaf, &      ! Distribute a state vector to model fields
         next_observation_pdaf, &      ! Provide time step of next observation
         prepoststep_ens_pdaf          ! User supplied pre/poststep routine

    ! Localization of state vector
    external :: init_n_domains_pdaf, & ! Provide number of local analysis domains
         init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
         g2l_state_pdaf, &             ! Get state on local analysis domain from global state
         l2g_state_pdaf                ! Update global state from state on local analysis domain

    ! Interface to PDAF-OMI for local and global filters
    external :: &
         init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
         obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
         init_dim_obs_l_pdafomi        ! Get dimension of obs. vector for local analysis domain


    ! *********************************
    ! *** Call assimilation routine ***
    ! *********************************

    ! Check  whether the filter is domain-localized
    call PDAF_get_localfilter(localfilter)

    if (localfilter == 1) then
       call PDAFomi_assimilate_local(collect_state_pdaf, &
            distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
            prepoststep_ens_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
            init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
            next_observation_pdaf, status_pdaf)
    elseif (localfilter == 0) then
       ! All global filters, except LEnKF
       call PDAFomi_assimilate_global(collect_state_pdaf, &
            distribute_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
            prepoststep_ens_pdaf, &
            next_observation_pdaf, status_pdaf)
    end if

! *** Query whether analysis step was performed
! *** This is also used to trigger the Euler time step for nemo_coupling='odir'
    call PDAF_get_assim_flag(assim_flag)

    ! Output into NEMO's ocean.output file
    if(assim_flag==1 .and. lwp) then
       write(numout,*) 
       write(numout,*) 'assimilate_pdaf : PDAF data assimilation was applied at step ', kt
       write(numout,*) '~~~~~~~~~~~'
    endif

    if (assim_flag==1) call MPI_Barrier(COMM_ensemble, MPIerr)

    if (coupling_nemo/='odir') assim_flag=0

    ! Check for errors during execution of PDAF
    if (status_pdaf /= 0) then
       write (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in PDAF_put_state - stopping! (PE ', mype_ens, ')'
       call abort_parallel()
    end if

  end subroutine assimilate_pdaf

end module assimilation_pdaf
