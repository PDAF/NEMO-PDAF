!> Collecting the state vector variables fro model fields
!!
!! The routine has to initialize the state vector of PDAF
!! from the fields of the model.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! **Calling Sequence**
!!
!!  - Called from:* `PDAFomi_assimilate_local`/`assimilation_pdaf` (as U_coll_state)
!!
subroutine collect_state_pdaf(dim_p, state_p)

  use mod_kind_pdaf

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector


  ! ********************
  ! Collect state vector
  ! ********************

  ! Nothing to be done for offline mode

end subroutine collect_state_pdaf
