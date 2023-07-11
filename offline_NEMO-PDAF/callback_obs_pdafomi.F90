!> PDAFOMI interface routines
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!! When adding an observation type, one has to add one module
!! obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!! In addition one has to add a call to the different routines include
!! in this file. It is recommended to keep the order of the calls
!! consistent over all files. 
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
subroutine init_dim_obs_pdafomi(step, dim_obs)

  use obs_ssh_mgrid_pdafomi, &
       only: assim_ssh_mgrid, init_dim_obs_ssh_mgrid
  use obs_sst_cmems_pdafomi, &
       only: assim_sst_cmems, init_dim_obs_sst_cmems

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step     !< Current time step
  integer, intent(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  integer :: dim_obs_ssh_mgrid         ! Observation dimension
  integer :: dim_obs_sst_cmems         ! Observation dimension


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_ssh_mgrid = 0
  dim_obs_sst_cmems = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called

  if (assim_ssh_mgrid) call init_dim_obs_ssh_mgrid(step, dim_obs_ssh_mgrid)
  if (assim_sst_cmems) call init_dim_obs_sst_cmems(step, dim_obs_sst_cmems)

  dim_obs = dim_obs_ssh_mgrid + dim_obs_sst_cmems

end subroutine init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
subroutine obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  use mod_kind_pdaf
  use obs_ssh_mgrid_pdafomi, &
       only: obs_op_ssh_mgrid
  use obs_sst_cmems_pdafomi, &
       only: obs_op_sst_cmems

  implicit none

! *** Arguments ***
  integer, intent(in) :: step                 !< Current time step
  integer, intent(in) :: dim_p                !< PE-local state dimension
  integer, intent(in) :: dim_obs              !< Dimension of full observed state
  real(pwp), intent(in)    :: state_p(dim_p)  !< PE-local model state
  real(pwp), intent(inout) :: ostate(dim_obs) !< PE-local full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of these calls is not relevant as the setup
  ! of the overall observation vector is defined by the
  ! order of the calls in init_dim_obs_pdafomi

   call obs_op_ssh_mgrid(dim_p, dim_obs, state_p, ostate)
   call obs_op_sst_cmems(dim_p, dim_obs, state_p, ostate)

end subroutine obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
subroutine init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  use mod_kind_pdaf
  use nemo_pdaf, &
       only: nwet
  use obs_ssh_mgrid_pdafomi, &
       only: assim_ssh_mgrid, init_dim_obs_l_ssh_mgrid
  use obs_sst_cmems_pdafomi, &
       only: assim_sst_cmems, init_dim_obs_l_sst_cmems

  implicit none

! *** Arguments ***
  integer, intent(in)  :: domain_p   !< Index of current local analysis domain
  integer, intent(in)  :: step       !< Current time step
  integer, intent(in)  :: dim_obs    !< Full dimension of observation vector
  integer, intent(out) :: dim_obs_l  !< Local dimension of observation vector

! *** Local variables ***

   if (nwet>0) then
      if (assim_ssh_mgrid) call init_dim_obs_l_ssh_mgrid(domain_p, step, dim_obs, dim_obs_l)
      if (assim_sst_cmems) call init_dim_obs_l_sst_cmems(domain_p, step, dim_obs, dim_obs_l)
   else
      dim_obs_l = 0
   end if

end subroutine init_dim_obs_l_pdafomi
