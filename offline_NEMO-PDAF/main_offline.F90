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
       only: dim_state, dim_state_p, assimilate_pdaf
  use initialize_offline, &
       only: initialize
  use initialize_pdaf, &      
       only: init_pdaf
  use utils_pdaf, &
       only: finalize_pdaf

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

! *** Initial Screen output ***
  initscreen: if (mype_world == 0) then

     write (*, '(/8x, a/)') '+++++ PDAF offline mode +++++'
     write (*, '(9x, a)') 'Data assimilation with PDAF'

     if (npes_world > 1) then
        write (*, '(21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     else
        write (*, '(21x, a/)') 'Running on 1 PE'
     end if

  end if initscreen

  
! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***

  call init_parallel_pdaf(comm_model)


! *** Initialize model information ***

  call initialize()


! *******************************
! ***      ASSIMILATION       ***
! *******************************

! *** Initialize PDAF ***

  call init_pdaf()

! *** Perform analysis ***

  call assimilate_pdaf()


! ********************
! *** Finishing up ***
! ********************

  call finalize_pdaf()

! *** Final screen output ***
  if (mype_world == 0) then
     ! *** Show allocated memory ***
     write (*, '(/22x, a)') 'Model - Memory overview'
     write (*, '(14x, 45a)') ('-', i=1, 45)
     write (*, '(25x, a)') 'Allocated memory  (MB)'
     write (*, '(17x, a, f12.3, a)') &
          'Model fields:', memcount_get(1, 'M'), ' MB'
     write (*, '(12x, a, f12.3, a)') &
          'forecast ensemble:', memcount_get(2, 'M'), ' MB'
  end if


! *** Terminate MPI
  call finalize_parallel()

end program MAIN_OFFLINE
