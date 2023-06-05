!> Main program for generiating covariance matrix
!!
!! This program is used to geenrate a covariance matrix
!! file for NEMO that can later be used to initialize an
!! ensemble for assimilation with PDAF.
!!
!! History:
!! - 2022-05 Lars Nerger, initial version based on PDAF offline code
!!
program MAIN_OFFLINE

  use parallel_pdaf, &        ! Parallelization
       only: MPI_COMM_WORLD, MPIerr, npes_world, mype_world, &
       init_parallel, finalize_parallel, init_parallel_pdaf, &
       comm_model
  use mod_memcount_pdaf, &    ! Memory counting
       only: memcount_ini, memcount_get
  use timer, &                ! Timings
       only: timeit, time_tot
  use assimilation_pdaf, &    ! Dimensions
       only: dim_state, dim_state_p
  use nemo_pdaf, &            ! NEMO-related variables and functions
       only: set_nemo_grid
    use statevector_pdaf, &   ! State vector functions
         only: setup_statevector

  implicit none

! Local variables
  integer :: i                 ! Counter


! **********************
! *** Initialize MPI ***
! **********************

  call init_parallel() ! initializes MPI


! ********************************
! ***      INITIALIZATION      ***
! ********************************

  ! Initialize memory counting and timers
  call memcount_ini(4)
  call timeit(3, 'ini')
  call timeit(1,'new')


! *** Initial Screen output ***
  initscreen: if (mype_world == 0) then

     write (*, '(/8x, a/)') '+++++ PDAF Generate Covar +++++'

     if (npes_world > 1) then
        write (*, '(16x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     else
        write (*, '(16x, a/)') 'Running on 1 PE'
     end if

  end if initscreen

  
! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***

  call init_parallel_pdaf(comm_model)


! ************************************************
! *** Specify state vector and state dimension ***
! ************************************************

  call timeit(2,'new')

  ! Initialize dimension information for NEMO grid
  call init_grid_dims()

  ! Initialize NEMO grid
  call set_nemo_grid()

  ! Setup state vector
  call setup_statevector(dim_state, dim_state_p)

  call timeit(2,'old')


! ************************************
! ***  GENERATE COVARIANCE MATRIX  ***
! ************************************

  call timeit(3,'new')
  call eofcovar()
  call timeit(3,'old')


! ********************
! *** Finishing up ***
! ********************

  call timeit(1,'old')

! *** Final screen output ***
  if (mype_world == 0) then
     ! *** Show allocated memory ***
     write (*, '(/22x, a)') 'Model - Memory overview'
     write (*, '(14x, 45a)') ('-', i=1, 45)
     write (*, '(25x, a)') 'Allocated memory  (MB)'
     write (*, '(17x, a, f12.3, a)') &
          'Model fields:', memcount_get(1, 'M'), ' MB'
     write (*, '(21x, a, f12.3, a)') &
          'eofcovar:', memcount_get(4, 'M'), ' MB'

     write (*, '(/17x, a)') 'Offline - Timing information'
     write (*, '(10x, 45a)') ('-', i=1, 45)
     ! Timing summary for assimilation
     write (*, '(19x, a, F11.3, 1x, a)') 'initialize model:', time_tot(2), 's'
     write (*, '(19x, a, F11.3, 1x, a)') 'generate covar:  ', time_tot(3), 's'
     write (*, '(19x, a, F11.3, 1x, a)') 'total run time:', time_tot(1), 's'

     write (*, '(/1x, a)') 'PDAF offline mode: END'
  end if


! *** Terminate MPI
  call finalize_parallel()

end program MAIN_OFFLINE
