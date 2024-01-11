!> Initialize PDAF
!!
!! This module contains the initialization routine for PDAF
!! `init_pdaf`. Here the ensemble is initialized and distributed
!! and the statevector and state variable information are computed.
!! 
!! Contributors:
!! The coupling and PDAF user codes for NEMO contains code parts
!! from different contributors. Namely
!! - Wibke Duesterhoeft-Wriggers, BSH, Germany
!! - Nicholas Byrne, NCEO and University of Reading, UK
!! - Yumeng Chen, NCEO and University of Reading, UK
!! - Yuchen Sun, AWI, Germany
!! - Lars Nerger, AWI, Germany
!!
module initialize_pdaf

   implicit none

contains

  !> This routine collects the initialization of variables for PDAF.
  !!
  !! The initialization routine `PDAF_init` is called
  !! such that the internal initialization of PDAF is performed.
  !! The initialization is used to set-up local domain and filter options
  !! such as the filter type, inflation, and localization radius.
  !!
  !! This variant is for the offline mode of PDAF.
  !!
  !! Much of the initialisation is read from a PDAF-specific namelist.
  !! This is performed in `read_config_pdaf`.
  !!
  subroutine init_pdaf()

    use mod_kind_pdaf
    use parallel_pdaf, &
         only: n_modeltasks, task_id, COMM_model, COMM_filter, &
         COMM_couple, COMM_ensemble, mype_ens, filterpe, abort_parallel
    use assimilation_pdaf, &
         only: dim_state, dim_state_p, screen, step_null, filtertype, &
         subtype, dim_ens, incremental, type_forget, forget, &
         type_trans, type_sqrt, delt_obs, locweight, type_ens_init, &
         type_central_state, type_hyb, hyb_gamma, hyb_kappa
    use nemo_pdaf, &
         only: set_nemo_grid
    use statevector_pdaf, &
         only: setup_statevector
    use utils_pdaf, &
         only: init_info_pdaf, read_config_pdaf
    use obs_sst_cmems_pdafomi, &
         only: assim_sst_cmems, rms_obs_sst_cmems, &
         lradius_sst_cmems, sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems
    use obs_ssh_mgrid_pdafomi, &
         only: assim_ssh_mgrid, rms_ssh_mgrid, &
         lradius_ssh_mgrid, sradius_ssh_mgrid
    use timer, only: timeit, time_temp

    implicit none

! *** Local variables ***
    integer :: i                 ! Counter
    integer :: filter_param_i(7) ! Integer parameter array for filter
    real    :: filter_param_r(3) ! Real parameter array for filter
    integer :: status_pdaf       ! PDAF status flag
    integer :: doexit, steps     ! Not used in this implementation
    real(pwp) :: timenow         ! Not used in this implementation

! *** External subroutines ***      
    external :: init_ens_pdaf    ! Ensemble initialization


! ***************************
! ***   Initialize PDAF   ***
! ***************************

    call timeit(2,'old')
    call timeit(3,'new')

    if (mype_ens == 0) then
       write (*, '(/a,1x,a)') 'NEMO-PDAF', 'INITIALIZE PDAF'
    end if


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

    ! *** IO options ***
    screen = 2  ! Write screen output (1) for output, (2) add timings

    ! *** Filter specific variables
    filtertype = 7    ! Type of filter
      !   (1) SEIK
      !   (2) EnKF
      !   (3) LSEIK
      !   (4) ETKF
      !   (5) LETKF
      !   (6) ESTKF
      !   (7) LESTKF
      !   (9) NETF
      !  (10) LNETF
      !  (11) LKNETF
      !  (12) PF
    dim_ens = n_modeltasks  ! Size of ensemble for all ensemble filters
      !   We use n_modeltasks here, initialized in init_parallel_pdaf
    subtype = 5       ! (5) Offline mode
    type_trans = 0    ! Type of ensemble transformation
      !   SEIK/LSEIK and ESTKF/LESTKF:
      !     (0) use deterministic omega
      !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
      !     (2) use product of (0) with random orthonormal matrix with
      !         eigenvector (1,...,1)^T
      !   ETKF/LETKF:
      !     (0) use deterministic symmetric transformation
      !     (2) use product of (0) with random orthonormal matrix with
      !         eigenvector (1,...,1)^T
    type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
      !   (0) fixed
      !   (1) global adaptive
      !   (2) local adaptive for LSEIK/LETKF/LESTKF
    forget = 1.0     ! Forgetting factor
    type_sqrt = 0     ! Type of transform matrix square-root
      !   (0) symmetric square root, (1) Cholesky decomposition
    type_hyb = 0      ! Type of hybrid weight:
    hyb_gamma =  1.0  ! Hybrid filter weight for state (1.0: LETKF, 0.0: LNETF)
    hyb_kappa = 30.0  ! Hybrid norm for using skewness and kurtosis


    ! ********************************************************************
    ! ***   Settings for ensemble init  - used in call-back routines   ***
    ! ********************************************************************

    type_ens_init = 2         ! Type of ensemble initialization
       !    (0) read snapshots from a single model file
       !    (1) read states from single ensemble file
       !    (2) read snapshots from separate model files
       !    (3) initialize ensemble from covariance matrix file
       !    (4) ensemble restart using fields from NEMO restart files

    type_central_state = 1    ! Type of central state of ensemble
       !    (0) mean of model snapshots
       !    (1) read from file
       !    (2) from NEMO field on model task using collect_state_pdaf 


    ! *********************************************************************
    ! ***   Settings for analysis steps  - used in call-back routines   ***
    ! *********************************************************************

    ! Which observations to assimilate
    assim_sst_cmems = .false.        ! Whether to assimilate SST data from CMEMS
    assim_ssh_mgrid = .false.        ! Whether to assimilate SSH data on model grid

    ! *** Localization settings
    locweight = 0     ! Type of localizating weighting
      !   (0) constant weight of 1
      !   (1) exponentially decreasing with SRADIUS
      !   (2) use 5th-order polynomial
      !   (3) regulated localization of R with mean error variance
      !   (4) regulated localization of R with single-point error variance

    ! Settings for SSH data on model grid
    rms_ssh_mgrid     = 0.8_8     ! Observation error stddev for SSH data on mopdel grid
    lradius_ssh_mgrid = 10000.0_8 ! Radius in km for lon/lat (or in grid points)
    sradius_ssh_mgrid = lradius_ssh_mgrid  ! Support radius for 5th-order polynomial
                                  ! or distance for 1/e for exponential weighting

    ! Settings for CMEMS satellite SST
    rms_obs_sst_cmems = 0.8_8     ! Observation error stddev for SST data from CMEMS
    mode_sst_cmems = 1            ! Observation mode for SST_CMEMS: 
                                  !  (0) linear interpolation to observation grid
                                  !  (1) super-obbing: average 4 observation values
    lradius_sst_cmems = 10000.0_8 ! Radius in km for lon/lat (or in grid points)
    sradius_sst_cmems = lradius_sst_cmems  ! Support radius for 5th-order polynomial
                                  ! or distance for 1/e for exponential weighting


! ******************************************
! *** Namelist reading and screen output ***
! ******************************************

    ! Read namelist file for PDAF
    call read_config_pdaf()

    ! Screen output for PDAF parameters
    if (mype_ens == 0) call init_info_pdaf()


! ************************************************
! *** Specify state vector and state dimension ***
! ************************************************

    ! Initialize dimension information for NEMO grid
    call set_nemo_grid(screen)

    ! Setup state vector
    call setup_statevector(dim_state, dim_state_p)


! *****************************************************
! *** Set communicator within which PDAF operates.  ***
! *****************************************************

    call PDAF_set_comm_pdaf(COMM_ensemble)


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! *****************************************************

    whichinit: if (filtertype /= 11) then
       ! *** All filters except LKNETF/EnKF/LEnKF ***
       filter_param_i(1) = dim_state_p ! State dimension
       filter_param_i(2) = dim_ens     ! Size of ensemble
       filter_param_i(3) = 0           ! Smoother lag (not implemented here)
       filter_param_i(4) = incremental ! Whether to perform incremental analysis
       filter_param_i(5) = type_forget ! Type of forgetting factor
       filter_param_i(6) = type_trans  ! Type of ensemble transformation
       filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
       filter_param_r(1) = forget      ! Forgetting factor

       call PDAF_init(filtertype, subtype, step_null, &
            filter_param_i, 7, &
            filter_param_r, 2, &
            COMM_model, COMM_filter, COMM_couple, &
            task_id, 1, filterpe, init_ens_pdaf, &
            screen, status_pdaf)
    else
       ! *** LKNETF ***
       filter_param_i(1) = dim_state_p ! State dimension
       filter_param_i(2) = dim_ens     ! Size of ensemble
       filter_param_i(3) = 0           ! Size of lag in smoother
       filter_param_i(4) = 0           ! Not used for NETF (Whether to perform incremental analysis)
       filter_param_i(5) = type_forget ! Type of forgetting factor
       filter_param_i(6) = type_trans  ! Type of ensemble transformation
       filter_param_i(7) = type_hyb    ! Type of hybrid weight
       filter_param_r(1) = forget      ! Forgetting factor
       filter_param_r(2) = hyb_gamma   ! Hybrid filter weight for state
       filter_param_r(3) = hyb_kappa   ! Normalization factor for hybrid weight 
     
       call PDAF_init(filtertype, subtype, step_null, &
            filter_param_i, 7, &
            filter_param_r, 3, &
            COMM_model, COMM_filter, COMM_couple, &
            task_id, 1, filterpe, init_ens_pdaf, &
            screen, status_pdaf)
    end if whichinit

    ! *** Check whether initialization of PDAF was successful ***
    if (status_pdaf /= 0) then
       write (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in initialization of PDAF - stopping! (PE ', mype_ens, ')'
       call abort_parallel()
    end if

    call timeit(3,'old')
    call timeit(4,'new')

  end subroutine init_pdaf

end module initialize_pdaf
