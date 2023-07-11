!> Define real precision
!!
!! This module defines the kind of real and length of character strings
!! for the PDAF call-back routines and interfaces. It is based on the NEMO
!! module `par_kind.F90`.
!!
MODULE mod_kind_pdaf

  IMPLICIT NONE
  SAVE

  !> double precision
  INTEGER, PARAMETER :: pdp = SELECTED_REAL_KIND(12, 307)
  !> double precision
  INTEGER, PARAMETER :: pwp = pdp
  !> standard string length
  INTEGER, PARAMETER :: lc = 256

END MODULE mod_kind_pdaf
