!> Setup parallelisation
!!
!! This module provides variables for the MPI parallelization
!! to be shared between model-related routines. There are variables
!! that are used in the model even without PDAF, and additional variables
!! that are only used if data assimilaion with PDAF is performed.
!! The initialization of communicators for execution with PDAF is
!! performed in `init_parallel_pdaf`.
!! 
module parallel_pdaf

  use mod_kind_pdaf

  implicit none
  save

  include 'mpif.h'

  ! Basic variables for model state integrations
  integer :: COMM_model         ! MPI communicator for model tasks
  integer :: mype_model         ! Rank in COMM_model
  integer :: npes_model         ! Size of COMM_model

  integer :: COMM_ensemble      ! Communicator for entire ensemble
  integer :: mype_ens           ! Rank in COMM_ensemble
  integer :: npes_ens           ! Size of COMM_ensemble

  ! Additional variables for use with PDAF
  integer :: n_modeltasks = 1   ! Number of parallel model tasks

  integer :: COMM_filter        ! MPI communicator for filter PEs 
  integer :: mype_filter        ! Rank in COMM_filter
  integer :: npes_filter        ! Size of COMM_filter

  integer :: COMM_couple        ! MPI communicator for coupling filter and model
  integer :: mype_couple        ! Rank in COMM_couple
  integer :: npes_couple        ! Size in COMM_couple

  integer :: mype_world         ! Rank in MPI_COMM_WORLD
  integer :: npes_world         ! Size in MPI_COMM_WORLD

  logical :: modelpe            ! Whether we are on a PE in a COMM_model
  logical :: filterpe           ! Whether we are on a PE in a COMM_filter
  integer :: task_id            ! Index of my model task (1,...,n_modeltasks)
  character(len=10) :: task_str ! Task ID as string
  integer :: MPIerr             ! Error flag for MPI
  integer :: MPIstatus(MPI_STATUS_SIZE)       ! Status array for MPI
  integer, allocatable :: local_npes_model(:) ! # PEs per ensemble

contains

!-------------------------------------------------------------------------------
!> Initialize MPI communicators for PDAF
!!
!! Split the MPI communicator initialised by XIOS into MODEL,
!! FILTER and COUPLE communicators, return MODEL communicator.
!!
   subroutine init_parallel_pdaf(mpi_comm)

#ifndef PDAF_OFFLINE
      use in_out_manager, only: cxios_context
#endif
      use timer, only: timeit

      !> Communicator after XIOS splitting
      integer, intent(inout) :: mpi_comm

      integer :: i, j                   !> Counters
      integer :: pe_index               !> Index of PE
      integer :: my_color, color_couple !> Variables for communicator-splitting
      integer :: dim_ens                !> Ensemble size / Number of model tasks
      integer :: tasks                  !> Number of taks for communicator splitting
      character(lc) :: nmlfile          !> Namelist file
      integer :: screen=1               !> Verbosity flag

      ! Number of ensemble members, supplied by PDAF namelist
      namelist /ensemble_nml/ dim_ens, screen

      call timeit(5,'ini')
      call timeit(5,'new')
      call timeit(1,'new')

      ! Read namelist for number of model tasks
      nmlfile = 'namelist_cfg.pdaf'

      open (20, file=nmlfile)
      read (20, NML=ensemble_nml)
      close (20)

      ! Online in case of online mode: dim_ens = number of parallel model tasks
      n_modeltasks = dim_ens
      tasks = dim_ens
#ifdef PDAF_OFFLINE
      ! offline mode always has only one model task
      tasks = 1
#endif

      ! ***              COMM_ENSEMBLE                ***
      ! *** Generate communicator for ensemble runs   ***
      ! *** only used to generate model communicators ***

      COMM_ensemble = mpi_comm
      call MPI_Comm_Size(COMM_ensemble, npes_ens, MPIerr)
      call MPI_Comm_Rank(COMM_ensemble, mype_ens, MPIerr)

      ! Initialize communicators for ensemble evaluations
      if (mype_ens == 0) then
         write (*, '(/a, 2x, a)') 'PDAF', 'Initialize communicators for assimilation with PDAF'
      end if

      ! Store # PEs per ensemble member. Used for info on PE 0 and for
      ! generation of model communicators on other PEs
      allocate (local_npes_model(tasks))
      local_npes_model = floor(real(npes_ens)/real(tasks))

      do i = 1, (npes_ens - tasks*local_npes_model(1))
         local_npes_model(i) = local_npes_model(i) + 1
      end do

      ! ***              COMM_MODEL               ***
      ! *** Generate communicators for model runs ***

      pe_index = 0
      doens1: do i = 1, tasks
         do j = 1, local_npes_model(i)
            if (mype_ens == pe_index) then
               task_id = i
               exit doens1
            end if
            pe_index = pe_index + 1
         end do
      end do doens1

      call MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
                          COMM_model, MPIerr)

      ! Re-initialize PE information according to model communicator
      call MPI_Comm_Size(COMM_model, npes_model, MPIerr)
      call MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

      if (screen > 1) then
         write (*, *) 'PDAF-MODEL: mype(w)= ', mype_ens, '; model task: ', task_id, &
            '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
      end if

      ! Init flag FILTERPE (all PEs of model task 1)
      if (task_id == 1) then
         filterpe = .true.
      else
         filterpe = .false.
      end if

      ! ***         COMM_FILTER                 ***
      ! *** Generate communicator for filter    ***

      if (filterpe) then
         my_color = task_id
      else
         my_color = MPI_UNDEFINED
      end if

      call MPI_Comm_split(COMM_ensemble, my_color, mype_ens, &
                          COMM_filter, MPIerr)

      ! Initialize PE information according to filter communicator
      if (filterpe) then
         call MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
         call MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)
      end if

      ! ***              COMM_COUPLE                 ***
      ! *** Generate communicators for communication ***
      ! *** between model and filter PEs             ***

      color_couple = mype_model + 1

      call MPI_Comm_split(COMM_ensemble, color_couple, mype_ens, &
                          COMM_couple, MPIerr)

      ! Initialize PE information according to coupling communicator
      call MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
      call MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

      if (screen > 0) then
         if (mype_ens == 0) then
            write (*, '(/18x, a)') 'PE configuration:'
            write (*, '(a, 2x, a6, a9, a10, a14, a13, /a, 2x, a5, a9, a7, a7, a7, a7, a7, /a, 2x, a)') &
               'Pconf', 'world', 'filter', 'model', 'couple', 'filterPE', &
               'Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
               'Pconf', '----------------------------------------------------------'
         end if
         call MPI_Barrier(COMM_ensemble, MPIerr)
         if (task_id == 1) then
            write (*, '(a, 2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
               'Pconf', mype_ens, mype_filter, task_id, mype_model, color_couple, &
               mype_couple, filterpe
         end if
         if (task_id > 1) then
            write (*, '(a, 2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
               'Pconf', mype_ens, task_id, mype_model, color_couple, mype_couple, filterpe
         end if
         call MPI_Barrier(COMM_ensemble, MPIerr)

         if (mype_ens == 0) write (*, '(/a)') ''

      end if

      ! ****************************************************
      ! *** Re-initialize model equivalent to COMM_model ***
      ! ****************************************************

      mpi_comm = COMM_model


#ifndef PDAF_OFFLINE
      ! Adapt XIOS contexts for ensemble
      write(task_str,'(I3.3)') task_id
      cxios_context = trim(cxios_context)//'_'//trim(task_str)
#endif

      call timeit(1,'old')
      call timeit(2,'new')

   end subroutine init_parallel_pdaf

!-------------------------------------------------------------------------------
!> Initialize the MPI execution environment.
!!
  subroutine init_parallel()

    implicit none

    integer :: i
  
    call MPI_INIT(i);
    call MPI_Comm_Size(MPI_COMM_WORLD,npes_world,i)
    call MPI_Comm_Rank(MPI_COMM_WORLD,mype_world,i)

    ! Initialize model communicator, its size and the process rank
    ! Here the same as for MPI_COMM_WORLD
    comm_model = MPI_COMM_WORLD
    npes_model = npes_world
    mype_model = mype_world
   
  end subroutine init_parallel

!-------------------------------------------------------------------------------
!> Finalize the MPI execution environment.
!!
  subroutine finalize_parallel()

    implicit none
    
    call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
    call  MPI_Finalize(MPIerr)

  end subroutine finalize_parallel

!-------------------------------------------------------------------------------
!> Terminate the MPI execution environment.
!!
  subroutine abort_parallel()

    call MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)

  end subroutine abort_parallel

end module parallel_pdaf
