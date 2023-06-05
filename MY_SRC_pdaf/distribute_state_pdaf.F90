!> Distributing the state vector variables, computing the
!> statevector increments
!!
!! This routine initializes statevector increments.
!! The increments are added to the model during the 
!! NEMO timestepping routine using NEMO's ASM module
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! **Calling Sequence**
!!
!!  - Called from: `PDAFomi_assimilate_local` (as U_dist_state)
!!
subroutine distribute_state_pdaf(dim_p, state_p)

  use mod_kind_pdaf
  use parallel_pdaf, &
       only: mype_ens
  use transforms_pdaf, &
       only: transform_field_mv
  use asminc_pdaf, &
       only: update_bkginc_pdaf

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: verbose                         ! Control verbosity


! ************************************************
! Compute increment for state vector variables ***
! for use in ASMINC module of NEMO             ***
! ************************************************

  if (mype_ens==0) then
     verbose = 1
  else
     verbose = 0
  end if

  ! Apply field transformations
  call transform_field_mv(2, state_p, 21, verbose)

  ! Update increment arrays for ASM
  call update_bkginc_pdaf(dim_p, state_p, verbose)

end subroutine distribute_state_pdaf
