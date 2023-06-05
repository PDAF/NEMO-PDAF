!$Id$
!BOP
!
! !MODULE:
module parser

! !DESCRIPTION:
! This module provides routine to parse command line
! arguments of different types. This version is for 
! use with MPI parallelization.
! By default, this routine uses the intrinsics 
! 'get\_command\_count' and 'get\_command\_argument' 
! that are define by the Fortran 2003 standard.
! If a compiler does not support these functions, you
! can use '-DF77' as a definition for the preprocessor.
! In this case the Fortran77 standard 'iargc()' and
! 'getarg()' are used.
!
! The module provides a generic subroutine to parse
! variables of type INTEGER, REAL, or CHARACTER
! (with length up to 100) from the command line.
!
! Usage:                      \begin{verbatim}
! SUBROUTINE PARSE(char(len=32) handle, variable)
!   The string 'handle' determines the name of    
!   the parsed variable.                          
!   Example: handle='iters' parses a variable     
!            specified on the command line by     
!            '-iters value'
!                                                 
!    Usage:                                       
!    CALL PARSE(handle, int_variable)             
!         Parses a variable of type integer       
!         whose name is given by the string       
!         handle.                                 
!                                                 
!    CALL PARSE(handle, real_variable)            
!         Parses a variable of type real          
!         whose name is given by the string       
!         handle.                                 
!                                                 
!    CALL PARSE(handle, character_variable)       
!         Parses a string variable of maxmimal    
!         length of 100 characters whose name is  
!         given by the string handle.             
!                                                 
!    CALL PARSE(handle, logical_variable)         
!         Parses a variable of type logical       
!         whose name is given by the string       
!         handle. In the command line it has      
!         to be specified as 'T' or 'F'.          
!                               \end{verbatim}
!
! !REVISION HISTORY:
! 2003-02 - Stephan Frickenhaus, Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  implicit none
  save
  
  include 'mpif.h'

! !PUBLIC MEMBER FUNCTIONS:
  public :: parse
  character(len=32), public :: handle  ! handle for command line parser
!EOP

  private
  character(len=100) :: str1, str2 
  integer :: i   
  integer :: mype, MPIerr
!   INTEGER,EXTERNAL :: iargc


  ! *** define interface ***
  interface parse
    module procedure parse_int
    module procedure parse_real
    module procedure parse_string
    module procedure parse_logical
  end interface

contains
  subroutine parse_int(handle, intvalue)
! ******************************
! *** Parse an integer value ***
! ******************************

! *** subroutine arguments ***    
    character(len=32), intent(in) :: handle
    integer,intent(inout) :: intvalue

! *** local variables ***
    character(len=32) :: string
    integer :: parsed_int
    logical :: modified

! *** Initialization ***
    call MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // trim(handle)
    modified = .false.
    
! *** Parsing ***
#ifdef F77
    write (*,*) 'PARSE for F77!!!!!!!!!!!!!!!'
    if (iargc() > 0) then 
       do i = 1, iargc() - 1 
          call getarg(i, str1) 
          call getarg(i + 1, str2) 
#else
    if (command_argument_count() > 0) then 
       do i = 1, command_argument_count() - 1 
          call get_command_argument(i, str1)
          call get_command_argument(i+1, str2)
#endif
          if (str1 == trim(string)) then
             read(str2, *) parsed_int
             modified = .true.
          end if
       enddo
    endif

! *** Finalize ***
    if (modified) then
       intvalue = parsed_int
!        IF (mype == 0) WRITE (*, '(2x, a, a, a, i)') &
       if (mype == 0) write (*, '(2x, a, a, a, i10)') &
            'PARSER: ', trim(handle), '=', parsed_int
    end if
  end subroutine parse_int

  subroutine parse_real(handle, realvalue)
! **************************
! *** Parse a real value ***
! **************************

! *** function arguments ***    
    character(len=32), intent(in) :: handle
    real, intent(inout) :: realvalue

! *** local variables ***
    character(len=32) :: string
    real :: parsed_real
    logical :: modified

! *** Initialize ***
    call MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // trim(handle)
    modified = .false.

! *** Parsing ***
#ifdef F77
    if (iargc() > 0) then 
       do i = 1, iargc() - 1 
          call getarg(i, str1) 
          call getarg(i + 1, str2) 
#else
    if (command_argument_count() > 0) then 
       do i = 1, command_argument_count() - 1 
          call get_command_argument(i, str1)
          call get_command_argument(i+1, str2)
#endif
          if (str1 == trim(string)) then
             read(str2, *) parsed_real
             modified = .true.
          end if
       enddo
    endif

! *** Finalize ***
    if (modified) then
       realvalue = parsed_real
       if (mype == 0) write (*, '(2x, a, a, a, es12.4)') &
            'PARSER: ', trim(handle), '=', parsed_real
    end if
  end subroutine parse_real


  subroutine parse_string(handle, charvalue)
! **********************
! *** Parse a string ***
! **********************

! *** function arguments ***    
    character(len=32), intent(in) :: handle
    character(len=*), intent(inout) :: charvalue

! *** local variables ***
    character(len=100) :: string
    character(len=100) :: parsed_string
    logical :: modified

! *** Initialize ***
    call MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // trim(handle)
    modified = .false.
    
! *** Parsing ***
#ifdef F77
    if (iargc() > 0) then 
       do i = 1, iargc() - 1 
          call getarg(i, str1) 
          call getarg(i + 1, str2) 
#else
    if (command_argument_count() > 0) then 
       do i = 1, command_argument_count() - 1 
          call get_command_argument(i, str1)
          call get_command_argument(i+1, str2)
#endif
          if (str1 == trim(string)) then
             read(str2, *) parsed_string
             modified = .true.
          end if
       enddo
    endif

! *** Finalize ***
    if (modified) then
       charvalue = parsed_string
       if (mype == 0) write (*, '(2x, a, a, a, a)') &
           'PARSER: ', trim(handle), '= ', trim(parsed_string)
    end if

  end subroutine parse_string

  subroutine parse_logical(handle, logvalue)
! ******************************
! *** Parse an logical value ***
! ******************************

! *** subroutine arguments ***    
    character(len=32), intent(in) :: handle
    logical, intent(inout) :: logvalue

! *** local variables ***
    character(len=32) :: string
    logical :: parsed_log
    logical :: modified

! *** Initialization ***
    call MPI_Comm_Rank(MPI_COMM_WORLD, mype, MPIerr)

    string = '-' // trim(handle)
    modified = .false.
    
! *** Parsing ***
#ifdef F77
    if (iargc() > 0) then 
       do i = 1, iargc() - 1 
          call getarg(i, str1) 
          call getarg(i + 1, str2) 
#else
    if (command_argument_count() > 0) then 
       do i = 1, command_argument_count() - 1 
          call get_command_argument(i, str1)
          call get_command_argument(i+1, str2)
#endif
          if (str1 == trim(string)) then
             read(str2, *) parsed_log
             modified = .true.
          end if
       enddo
    endif

! *** Finalize ***
    if (modified) then
       logvalue = parsed_log
       if (mype == 0) write (*, '(2x, a, a, a, l1)') &
            'PARSER: ', trim(handle), '=', parsed_log
    end if
  end subroutine parse_logical

end module parser
