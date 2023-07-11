!> Set number of local analysis domains
!!
!! The routine is called in `PDAF_X_update`
!! at the beginning of the analysis step before
!! the loop through all local analysis domains.
!! It has to set the number of local analysis
!! domains for the PE-local domain.
!!
!! This code is for NEMO-PDAF
!!
!! - Called from: `PDAFomi_assimilate_local`/`assimilation_pdaf`
!
subroutine init_n_domains_pdaf(step, n_domains_p)

  use nemo_pdaf, &
       only: nwet

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step        !< Current time step
  integer, intent(out) :: n_domains_p !< PE-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************

  ! *******************************************
  !
  ! The number of local domains is defined as
  ! the number of grid points at the surface
  ! where tmask is 1 i.e. horizontal localization
  ! is used, and land points are ignored.
  !
  ! *******************************************

  n_domains_p = nwet

end subroutine init_n_domains_pdaf
