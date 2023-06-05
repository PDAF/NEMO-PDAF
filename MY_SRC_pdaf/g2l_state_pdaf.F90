!>Restrict a model state to a local analysis domain
!!
!! The routine is called during the loop over all
!! local analysis domains in `PDAF_X_update`
!! before the analysis on a single local analysis
!! domain. It has to initialize elements of the
!! state vector for the local analysis domains from
!! the PE-local full state vector.
!!
subroutine g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

  use mod_kind_pdaf
  use assimilation_pdaf, &
       only: id_lstate_in_pstate

  implicit none

! *** Arguments ***
  integer, intent(in) :: step              !< Current time step
  integer, intent(in) :: domain_p          !< Current local analysis domain
  integer, intent(in) :: dim_p             !< PE-local full state dimension
  integer, intent(in) :: dim_l             !< Local state dimension
  real(pwp), intent(in)  :: state_p(dim_p) !< PE-local full state vector 
  real(pwp), intent(out) :: state_l(dim_l) !< State vector on local analysis domain

! *** Local variables ***
  integer :: i                             ! Counter


  ! *************************************
  ! *** Initialize local state vector ***
  ! *************************************

  ! This generic initialization uses id_lstate_in_pstate
  ! which has been initialied in init_dim_l_pdaf

  do i = 1, dim_l
     state_l(i) = state_p(id_lstate_in_pstate(i))
  enddo

end subroutine g2l_state_pdaf
