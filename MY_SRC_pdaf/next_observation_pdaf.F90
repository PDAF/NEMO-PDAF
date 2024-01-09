!> Determining the Next Analysis Step
!!
!! The subroutine is called before each forecast phase
!! by `PDAF_get_state`. It has to initialize the number
!! of time steps until the next available observation
!! (`nsteps`). It indicates if the data assimilation process
!! is completed such that the ensemble loop in the model
!! routine can be exited.
!! 
!! The routine is called by all processes.
!! 
!!  - Called from: `init_pdaf/PDAF_get_state` (as U_next_obs)
!! 
subroutine next_observation_pdaf(stepnow, nsteps, doexit, time)

  use mod_kind_pdaf
  use assimilation_pdaf, &
       only: delt_obs
  use parallel_pdaf, &
       only: mype_ens
  use nemo_pdaf, &
       only: nitend, nit000
  use asminc_pdaf, &
       only: do_asmiau, steps_asmiau, &
       store_asm_step_pdaf, update_asm_step_pdaf
#if defined key_top
  use asminc_pdaf, &
       only: do_bgciau, steps_bgciau
#endif

  implicit none

! *** Arguments ***   
  integer, intent(in)  :: stepnow !< Number of the current time step
  integer, intent(out) :: nsteps  !< Number of time steps until next obs
  integer, intent(out) :: doexit  !< Whether to exit forecasting (1 for exit)
  real(pwp), intent(out) :: time  !< Current model (physical) time


! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************

  ! Not used in this implementation
  time = 0.0
  doexit = 0

  if (stepnow + delt_obs + nit000 - 1 <= nitend) then
     ! *** During the assimilation process ***

     ! The time step settings here account for the fact that for
     ! direct initialization one needs one time steps afterwards
     ! to appliy the increment, while for IAU one needs as many
     ! time steps as set for IAU. A DA run has to include these
     ! steps before writing restart files as otherwise the increment
     ! will not be fully applied. Here, we subtract this number
     ! of steps to ensure the application of the increment. An
     ! alternative would be to compute the first analysis step
     ! before the time stepping and not to perform an analysis
     ! step at the very end.
     
     if (stepnow == 0) then
        ! First analysis step 

        ! Apply IAU
#if defined key_top
        if (do_asmiau .or. do_bgciau) then
           if (steps_asmiau>= steps_bgciau) then
              nsteps = delt_obs - steps_asmiau
           else
              nsteps = delt_obs - steps_bgciau
           end if
#else
        if (do_asmiau) then
           nsteps = delt_obs - steps_asmiau
#endif
        else
           ! Direct initialization with increments
           nsteps = delt_obs-1     ! Analysis step one step before end of day
        endif
     else
        nsteps = delt_obs       ! Follow-up analysis steps daily
     end if

     if (mype_ens == 0) write (*, '(a, i7, 3x, a, i7)') &
          'NEMO-PDAF', stepnow, 'Next observation at time step', stepnow + nsteps

     ! Update analysis step information for NEMO-ASM
     if (stepnow == 0) then
        ! At initial time - apply increments after first analysis step
        call store_asm_step_pdaf(stepnow+nsteps, 0)
     else
        ! analysis step - apply increments after current analysis step
        call store_asm_step_pdaf(stepnow, 1)
     end if

  else
     ! *** End of assimilation process ***

     ! Set nsteps so stepnow+delt_bs>nitend to ensure that PDAF calls distribute_state
     nsteps = delt_obs        

     if (mype_ens == 0) write (*, '(a, i7, 3x, a)') &
          'NEMO-PDAF', stepnow, 'No more observations - end assimilation and complete forecast'

     ! Update analysis step information for NEMO-ASM
     call store_asm_step_pdaf(stepnow, 1)
  end if

end subroutine next_observation_pdaf
