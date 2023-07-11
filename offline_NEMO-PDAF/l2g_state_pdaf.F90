!> Initialize full state from local state vector
!!
!! The routine is called during the loop over all
!! local analysis domains in `PDAFomi_assimilate_local`
!! after the analysis and ensemble transformation
!! on a single local analysis domain. It has to
!! initialize elements of the PE-local full state
!! vector from the provided analysis state vector
!! on the local analysis domain.
!!
subroutine l2g_state_pdaf(step, domain_p_all, dim_l, state_l, dim_p, state_p)

  use mod_kind_pdaf
  use assimilation_pdaf, &
       only: id_lstate_in_pstate
  use statevector_pdaf, &
       only: n_fields, sfields, sfields_l
  use nemo_pdaf, &
       only: nwet

  implicit none

! *** Arguments ***
  integer, intent(in) :: step                !< Current time step
  integer, intent(in) :: domain_p_all        !< Current local analysis domain
  integer, intent(in) :: dim_l               !< Local state dimension
  integer, intent(in) :: dim_p               !< PE-local full state dimension
  real(pwp), intent(in)    :: state_l(dim_l) !< State vector on local analysis domain
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local full state vector 

! *** local variables *** 
  integer :: i, ifield           ! Counters


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  do ifield = 1, n_fields
     do i = sfields_l(ifield)%off+1, sfields_l(ifield)%off + sfields_l(ifield)%dim
        state_p(id_lstate_in_pstate(i)) = state_l(i)
     end do
  end do

end subroutine l2g_state_pdaf
