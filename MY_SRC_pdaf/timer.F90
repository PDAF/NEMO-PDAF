!$Id: timer.F90 266 2019-11-28 12:49:57Z lnerger $
!BOP
!
! !MODULE:
module timer

! !DESCRIPTION: 
! This module provides methods to perform timings of 
! parts of a program execution. It uses the intrinsic 
! function SYSTEM\_CLOCK.
!
! !REVISION HISTORY:
! 2000-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  implicit none
  save
  
  public :: timeit, time_tot, time_temp
!EOP

  private
  integer :: t_rate
  integer, allocatable :: t_start(:), t_end(:)
  real, allocatable    :: t_total(:), t_temp(:)

contains
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: timeit - Initialize Counters and time regions
!
! !INTERFACE: timeit()
  subroutine timeit(timerID, operation)

! !DESCRIPTION:
! Subroutine to initialize counters and to perform timing of a region
! specified by timerID.
! Usage:\\
!   CALL PDAF\_timeit(N,'ini') - Allocates and initializes N counters\\
!   CALL PDAF\_timeit(M,'new') - Start timing region for counter M\\
!   CALL PDAF\_timeit(M,'old') - End timing region for counter M\\
!   CALL PDAF\_timeit(M,'fin') - Finalized and deallocates all counters\\

! !USES:
    implicit none

! !ARGUMENTS:
    integer, intent(in) :: timerID             ! ID of timer
    character(len=3), intent(in) :: operation  ! Requested operation 
!EOP

    ! Initialize timers
    if (operation == 'ini') then
       if ( .not. (allocated(t_start))) then
          allocate(t_start(timerID), t_end(timerID))
          allocate(t_total(timerID), t_temp(timerID))
       end if
        
       t_total = 0.0
    end if
    
    ! Begin timing region
    if (operation == 'new') then
       call system_clock(t_start(timerID))
    end if

    ! End timing region
    if (operation == 'old') then
       call system_clock(t_end(timerID), t_rate)
       t_temp(timerID) = real(t_end(timerID) - t_start(timerID)) &
            / real(t_rate)
       t_total(timerID) = t_total(timerID) + real(t_end(timerID) - &
            t_start(timerID)) / real(t_rate)
    end if
    
    ! Finalize timers
    if (operation == 'fin') then
       deallocate(t_start, t_end)
       deallocate(t_total, t_temp)
    end if
    
  end subroutine timeit

!-------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: time_temp - Read out timers for last timing interval
!
! !INTERFACE: time_temp()
  real function time_temp(timerID)

! !DESCRIPTION:
! Read out the value of the timer in seconds for the last 
! passage of the timing region defined by timerID.

! !USES:
    implicit none

! !ARGUMENTS:
    integer, intent(in) :: timerID             ! ID of timer
!EOP

    time_temp = t_temp(timerID)

  end function time_temp

!-------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: PDAF_time_tot - Read out total time of a timing region.
!
! !INTERFACE: time_tot()
    real function time_tot(timerID)

! !DESCRIPTION:
! Read out the accumulated value of the timer in seconds
! for the timing region define by timerID.

! !USES:
    implicit none

! !ARGUMENTS:
    integer, intent(in) :: timerID             ! ID of timer
!EOP

    time_tot = t_total(timerID)

  end function time_tot

end module timer
