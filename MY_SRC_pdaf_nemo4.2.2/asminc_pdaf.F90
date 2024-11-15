!> Module for interfacing between PDAF and NEMO's ASMINC module
!!
!! This module provides the interfacing with NEMO's
!! ASMINC module. Settings of ASMINC are set on the
!! basis of DA settings for PDAF. both the direct
!! initialization and IAU of ASMINC can be used.
!!
module asminc_pdaf

  use mod_kind_pdaf
  use statevector_pdaf, &
       only: update_ssh, update_temp, update_salt, update_vel, &
       sfields, id
#if defined key_top_asm
  use statevector_pdaf, &
       only: n_trc, sv_trc
  use nemo_pdaf, &
       only: jptra
#endif
  use parallel_pdaf, &
       only: mype_ens
  use par_oce, &
       only: jpi, jpj, jpk, jp_tem, jp_sal
  use in_out_manager, &
       only: nit000, nitend
#if defined key_qco   ||   defined key_linssh
  use stpmlf, &
       only: Nbb, Nnn, Nrhs
#else
  use step, &
       only: Nbb, Nnn, Nrhs
#endif
  use asminc, &
       only: ln_bkgwri, ln_trainc, ln_dyninc, ln_sshinc, &
       ln_asmdin, ln_asmiau, nitbkg, nitdin, nitiaustr, &
       nitiaufin, niaufn, ln_salfix, salfixmin, nn_divdmp, &
       tra_asm_inc, dyn_asm_inc, ssh_asm_inc, wgtiau
  use asmpar, &
       only: nitdin_r, nitiaustr_r, nitiaufin_r
#if defined key_top_asm
  use asminc, &
       only: n_update_trc, ids_update_trc, bgc_bkginc, &
       nitdinbgc, nitibgcstr, nitibgcfin, niaufnbgc, &
       ln_trcinc, ln_bgcdin, ln_bgciau, wgtiaubgc, &
       nitdinbgc_r, nitibgcstr_r, nitibgcfin_r
#endif

  implicit none
  save
  public

  ! Variables that control the behavior of increments
  ! Next to these, the DA parameters set for PDAF influence the increments
  logical :: do_asmiau = .false.      ! Whether to apply IAU for ocean physics; if FALSE, DIN is used
  integer :: steps_asmiau = 1         ! Number of steps for physics IAU 
  integer :: shape_asmiau = 0         ! Shape of physics IAU function: (0) constant; (1) hat (NEMO: niaufn)
  real(pwp), allocatable :: iauweight(:)   ! Weights for IAU
  integer :: iter_divdmp = 0          ! Number of iterations of divergence damping operator

#if defined key_top
  logical :: do_bgciau = .false.      ! Whether to apply IAU for BGC; if FALSE, DIN is used
  integer :: steps_bgciau = 1         ! Number of steps for BGC IAU
  integer :: shape_bgciau = 0         ! Shape of BGC IAU function: (0) constant; (1) hat (NEMO: niaufnbgc)
  real(pwp), allocatable :: bgciauweight(:)   ! Weights for IAU for BGC 
#endif

! *** Local variables ***
  integer, private :: next_inc = 0

contains

!> Routine to initialise variables for NEMO's ASMINC module
!!
!! The routine initializes variables for the ASMINC module
!! of NEMO. Without PDAF they would be read from a namelist.
!! The routine also initializes the increment array for the
!! BGC data assimilation uypdate.
!!
  subroutine asm_inc_init_pdaf(delt_obs)

    implicit none

! *** Arguments ***
    integer, intent(in) :: delt_obs

! *** Local variables
    integer :: i             ! Counter
    real(pwp) :: totwgt      ! Accumulated IAU weight for checking
#if defined key_top_asm
    integer :: id_trc        ! Index of BGC variable in trc array
    integer :: id_var        ! Index of BGC variable in state vector
#endif
     


! *******************************************
! *** Initialize Settings for NEMO ASMINC ***
! *******************************************

    ! Set initial increment step
    next_inc = delt_obs

    ! Set switches and parameters for ASMINC - physics

    ln_bkgwri  = .false.   !  Logical switch for writing out background state
    if (.not. (update_ssh .or. update_temp .or. update_salt .or. update_vel)) then
       ln_trainc  = .false.   !  Logical switch for applying tracer increments
       ln_dyninc  = .false.   !  Logical switch for applying velocity increments
       ln_sshinc  = .false.   !  Logical switch for applying SSH increments
       ln_asmdin  = .false.      !  Logical switch for Direct Initialization (DI)
       ln_asmiau  = .false.   !  Logical switch for Incremental Analysis Updating (IAU)
    else
       if (update_temp .or. update_salt) then
          ln_trainc  = .true.    !  Logical switch for applying tracer increments
       else
          ln_trainc  = .false.   !  Logical switch for applying tracer increments
       end if
       if (update_vel) then
          ln_dyninc  = .true.    !  Logical switch for applying velocity increments
       else
          ln_dyninc  = .false.   !  Logical switch for applying velocity increments
       end if
       if (update_ssh) then
          ln_sshinc  = .true.    !  Logical switch for applying SSH increments
       else
          ln_sshinc  = .false.   !  Logical switch for applying SSH increments
       end if
       if (do_asmiau) then
          ln_asmiau  = .true.    !  Logical switch for Incremental Analysis Updating (IAU)
          ln_asmdin  = .false.   !  Logical switch for Direct Initialization (DI)
       else
          ln_asmiau  = .false.   !  Logical switch for Incremental Analysis Updating (IAU)
          ln_asmdin  = .true.    !  Logical switch for Direct Initialization (DI)
       end if
    end if

    ! The step settings are initialized here and later updated in cycled DA
    nitbkg    = 0                      !  Timestep of background in [0,nitend-nit000-1]
    nitdin    = delt_obs               !  Timestep of background for DI in [0,nitend-nit000-1]
    nitiaustr = delt_obs               !  Timestep of start of BGC IAU interval
    nitiaufin = delt_obs+steps_asmiau  !  Timestep of end of BGC IAU interval
    niaufn    = shape_asmiau           !  Type of IAU weighting function
    ln_salfix = .false.                !  Logical switch for ensuring that the sa > salfixmin
    salfixmin = -9999                  !  Minimum salinity after applying the increments
    nn_divdmp = iter_divdmp            !  Number of iterations of divergence damping operator


#if defined key_top_asm

    ! Switches and parameters for ASM BGC increments

    ln_trcinc  = .true.        !  Logical switch for applying BGC increments

    if (do_bgciau) then        ! Switch between IAU and direct insertion
       ln_bgciau  = .true.     !  Logical switch for Incremental Analysis Updating (IAU)
       ln_bgcdin = .false.     !  Logical switch for applying DI for BGC
    else
       ln_bgciau  = .false.    !  Logical switch for Incremental Analysis Updating (IAU)
       ln_bgcdin = .true.      !  Logical switch for applying DI for BGC
    end if
    niaufnbgc  = shape_bgciau  !  Type of BGC IAU weighting function
     
    ! The step settings are initialized here and later updated in cycled DA
    nitdinbgc  = delt_obs                !  Timestep for DI for BGC
    nitibgcstr = delt_obs                !  Timestep of start of BGC IAU interval
    nitibgcfin = delt_obs+steps_bgciau   !  Timestep of end of BGC IAU interval
#endif


! *****************************************
! *** Initialize weights vector for IAU ***
! *****************************************

    if (ln_asmiau) then
       allocate(iauweight(steps_asmiau))
       call init_iauweight(iauweight, steps_asmiau, shape_asmiau)
    end if

#if defined key_top_asm
    if (ln_bgciau) then
       allocate(bgciauweight(steps_bgciau))
       call init_iauweight(bgciauweight, steps_bgciau, shape_bgciau)
    end if
#endif


! *************************************************
! *** Prepare increment array for BCG variables ***
! *************************************************

#if defined key_top_asm

    ! Count number of prognostic BGC fields that are updated by the DA
    allocate(ids_update_trc(jptra))
    do id_trc = 1, jptra

       if (sv_trc(id_trc)) then

          id_var=id%trc(id_trc)

          if (sfields(id_var)%update) then
             n_update_trc = n_update_trc + 1
             ids_update_trc(n_update_trc) = id_trc
          end if
       end if
    end do

    ! Allocate increment array for BGC variables
    allocate(bgc_bkginc(jpi, jpj, jpk, n_update_trc))
    bgc_bkginc = 0.0_pwp

    ! Update ln_trcinc if no BGC updates are done by PDAF
    if (n_update_trc==0) ln_trcinc  = .false.    !  Logical switch for applying tracer increments

#endif


! *********************
! *** Screen output ***
! *********************

    if (mype_ens==0) then
       write (*,'(/a,4x,a)') 'NEMO-PDAF', '******** Setup for NEMO ASM ********'
       if (update_temp .or. update_salt .or. update_ssh .or. update_vel) then 
          write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Apply increment for NEMO physics'
          if (update_temp) write (*,'(a,10x,a)') 'NEMO-PDAF', '--- update temperature'
          if (update_salt) write (*,'(a,10x,a)') 'NEMO-PDAF', '--- update salinity'
          if (update_ssh)  write (*,'(a,10x,a)') 'NEMO-PDAF', '--- update SSH'
          if (update_vel)  write (*,'(a,10x,a)') 'NEMO-PDAF', '--- update velocities'
          if (nn_divdmp>0) &
               write (*,'(a,8x,a,i9)') 'NEMO-PDAF', '--- Apply divergence damping, iterations', nn_divdmp
          if (ln_asmdin)   write (*,'(a,10x,a)') 'NEMO-PDAF', '--- Use DIN for NEMO fields'
          if (ln_asmiau) then
             write (*,'(a,8x,a)') 'NEMO-PDAF', '--- Use IAU for NEMO fields'
             write (*,'(a,8x,a,i9)') 'NEMO-PDAF', '--- Number of IAU steps', steps_asmiau
             write(*,'(a,10x,a)') 'NEMO-PDAF', '--- IAU weights'
             write(*,'(a,18x,a,5x,a)') 'NEMO-PDAF', 'IAU step','IAU  weight'

             totwgt = 0.0
             do i = 1, steps_asmiau
                totwgt = totwgt + iauweight(i)
                write(*,'(a,14x,i6,7x,f12.6)') 'NEMO-PDAF', i, iauweight(i) 
             end do
             write(*,'(a,14x,a,f10.6)') 'NEMO-PDAF', 'Time-integrated weight = ', totwgt
          end if
       else
          write (*,'(a,4x,a)') 'NEMO-PDAF', '--- No increment for NEMO physics'
       end if

#if defined key_top_asm
       if (ln_trcinc) then
          write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Apply increment for BGC'
          write (*,'(a,4x,a,i5)') 'NEMO-PDAF', '--- Number of updated BCG fields:', n_update_trc
          if (ln_bgciau) then
             write (*,'(a,8x,a)') 'NEMO-PDAF', '--- Use IAU for BGC fields'
             write (*,'(a,8x,a,i9)') 'NEMO-PDAF', '--- Number of BGC IAU steps', steps_bgciau
             write(*,'(a,10x,a)') 'NEMO-PDAF', '--- IAU weights'
             write(*,'(a,18x,a,5x,a)') 'NEMO-PDAF', 'IAU step','IAU  weight'

             totwgt = 0.0
             do i = 1, steps_bgciau
                totwgt = totwgt + bgciauweight(i)
                write(*,'(a,14x,i6,7x,f12.6)') 'NEMO-PDAF', i, bgciauweight(i) 
             end do
             write(*,'(a,14x,a,f10.6)') 'NEMO-PDAF', 'Time-integrated weight = ', totwgt
          else
             write (*,'(a,8x,a)') 'NEMO-PDAF', '--- Use DIN for BGC fields'
          end if

       else
          write (*,'(a,4x,a,i5)') 'NEMO-PDAF', '--- No increment for BCG fields'
       end if
#endif
       write (*,'(a,4x,a)') 'NEMO-PDAF', '************************************'
    end if

    ! Update increment step
    call update_asm_step_pdaf(0)

  end subroutine asm_inc_init_pdaf



!-------------------------------------------------------------------------------
!> Routine to store analysis step for NEMO-ASM
!!
!! For direct initialization, the increment is applied on the
!! time step after the analysis update is computed by PDAF
  subroutine store_asm_step_pdaf(nextinc, initial)

    implicit none

! *** Arguments ***
    integer, intent(in) :: nextinc      ! Time of (starting) next assimilation increment
    integer, intent(in) :: initial      ! Whether the routine is called at inital time

    ! Store step of (start) of next increments
    next_inc = nextinc + 1

    call update_asm_step_pdaf(initial)

  end subroutine store_asm_step_pdaf



!-------------------------------------------------------------------------------
!> Routine to update analysis step values for NEMO-ASM
!!
!! This routine is called in assimilate_pdaf to
!! update the time step information for applying
!! assimilation increments.
!!
  subroutine update_asm_step_pdaf(initial)

    implicit none

! *** Arguments ***
    integer, intent(in) :: initial

! *** Local variables ***
    integer :: i, cnt            ! Counter


! *** Settings for ASMINC physics ***

    nitdin       = next_inc                           ! Time step of the background state for direct initialization
    nitdin_r     = nitdin    + nit000 - 1             ! Background time for DI referenced to nit000
    nitiaustr   = next_inc                            ! Timestep of start of IAU interval
    nitiaufin   = nitiaustr+steps_asmiau - 1          ! Timestep of end of IAU interval
    nitiaustr_r = nitiaustr + nit000 - 1              ! Start of IAU interval referenced to nit000
    nitiaufin_r = nitiaufin + nit000 - 1              ! End of IAU interval referenced to nit000

    if (mype_ens==0) then
       if (ln_trainc .or. ln_dyninc .or. ln_sshinc) then
          if (ln_asmiau) then
             write (*,'(a,5x,a,2i7)') 'NEMO-PDAF', '--- set IAU steps for ASMINC: ', nitiaustr_r, nitiaufin_r
          else
             write (*,'(a,5x,a,i7)') 'NEMO-PDAF', '--- set DIN step for ASMINC: ', nitdin_r
          end if
       end if
    end if

    ! Set weights in the full IAU weight array
    if (ln_asmiau .and. initial > 0) then

       if (allocated(wgtiau)) then
          if (mype_ens==0) &
               write (*,'(a,5x,a,2i7)') 'NEMO-PDAF', '--- update IAU weight array'

          cnt = 1
          wgtiau(:) = 0.0
          do i = nitiaustr, nitiaufin
             wgtiau(i) = iauweight(cnt)
             cnt = cnt + 1
          end do
       else
          write (*,*) 'PDAF-ERROR: wgtiau not allocated; is program compiled with key_asminc?'
       end if
    end if


#if defined key_top_asm
! *** Settings for ASMINC BGC ***

    nitdinbgc    = next_inc                           ! Time step of direct init for BGC
    nitdinbgc_r  = nitdinbgc    + nit000 - 1          ! Background time for DI referenced to nit000
    nitibgcstr   = next_inc                           ! Timestep of start of BGC IAU interval
    nitibgcfin   = nitibgcstr+steps_bgciau - 1        ! Timestep of end of BGC IAU interval
    nitibgcstr_r = nitibgcstr + nit000 - 1            ! Start of BGC IAU interval referenced to nit000
    nitibgcfin_r = nitibgcfin + nit000 - 1            ! End of BGC IAU interval referenced to nit000

    if (mype_ens==0) then
       if (ln_trcinc) then
          if (ln_bgciau) then
             write (*,'(a,5x,a,2i7)') 'NEMO-PDAF', '--- set BGC IAU steps for ASMINC: ', nitibgcstr_r, nitibgcfin_r
          else
             write (*,'(a,5x,a,i7)') 'NEMO-PDAF', '--- set BGC DIN step for ASMINC: ', nitdinbgc_r
          end if
       end if
    end if

    ! Set weights in the full BGC IAU weight array
    if (ln_bgciau .and. initial > 0) then

       if (allocated(wgtiaubgc)) then
          if (mype_ens==0) &
               write (*,'(a,5x,a,2i7)') 'NEMO-PDAF', '--- update BGC IAU weight array'

          cnt = 1
          wgtiaubgc(:) = 0.0
          do i = nitibgcstr, nitibgcfin
             wgtiaubgc(i) = bgciauweight(cnt)
             cnt = cnt + 1
          end do
       else
          write (*,*) 'PDAF-ERROR: wgtiaubgc not allocated; is program compiled with key_asminc?'
       end if
    end if
#endif

  end subroutine update_asm_step_pdaf



!-------------------------------------------------------------------------------
!> Routine to update the bkginc arrays of NEMO-ASM
!!
!! This routine updates the increment arrays for the
!! ASMINC module of NEMO.
!!
  subroutine update_bkginc_pdaf(dim_p, state_p, verbose)

    use nemo_pdaf, &
         only: ni_p, nj_p, nk_p, i0, j0
    use transforms_pdaf, &
         only: state2field_inc
    use oce, &
         only: ssh, ts, uu, vv
    use dom_oce, &
         only: tmask, umask, vmask
    use lbclnk, &
         only: lbc_lnk
    use asminc, &
         only: ssh_bkg, t_bkg, s_bkg, u_bkg, v_bkg, &
         ssh_bkginc, t_bkginc, s_bkginc, u_bkginc, v_bkginc
#if defined key_top_asm
    use trc, &
         only: tr
#endif

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p               !< PE-local state dimension
    real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector
    integer, intent(in) :: verbose             !< Verbosity flag

! *** Local variables ***
    integer :: i_bgcinc, i_tr     ! Indices
    integer :: id_var              ! Index


! *************************************************
! *** Prepare increment arrays for NEMO physics ***
! *************************************************

    if (verbose==1) write (*,'(a,4x,a)') 'NEMO-PDAF', 'Distribute state increment'

    physics: if (update_temp .or. update_salt .or. update_ssh .or. update_vel) then

       ! Ensure that the increment arrays are allocated and set to zero
       if (.not. allocated(ssh_bkginc)) allocate(ssh_bkginc(jpi,jpj))
       if (.not. allocated(t_bkginc)) allocate(t_bkginc(jpi,jpj,jpk))
       if (.not. allocated(s_bkginc)) allocate(s_bkginc(jpi,jpj,jpk))
       if (.not. allocated(u_bkginc)) allocate(u_bkginc(jpi,jpj,jpk))
       if (.not. allocated(v_bkginc)) allocate(v_bkginc(jpi,jpj,jpk))
       ssh_bkginc = 0.0_pwp
       t_bkginc = 0.0_pwp
       s_bkginc = 0.0_pwp
       u_bkginc = 0.0_pwp
       v_bkginc = 0.0_pwp

       ! SSH
       if (id%ssh > 0 .and. update_ssh) then
          call state2field_inc(state_p, ssh(1+i0:ni_p+i0, 1+j0:nj_p+j0, Nnn), &
               ssh_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0), sfields(id%ssh)%off, sfields(id%ssh)%ndims)

          ssh_bkginc(:,:) = ssh_bkginc(:,:) * tmask(:,:,1)

          ! Fill halo regions
          call lbc_lnk('distribute_state_pdaf', ssh_bkginc, 'T', 1.)
       end if

       ! T
       if (id%temp > 0 .and. update_temp) then
          call state2field_inc(state_p, &
               ts(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem, Nnn), &
               t_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
               sfields(id%temp)%off, sfields(id%temp)%ndims)
       end if
       t_bkginc(:,:,:) = t_bkginc(:,:,:) * tmask(:,:,:)

       ! S
       if (id%salt > 0 .and. update_salt) then
          call state2field_inc(state_p, &
               ts(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_sal, Nnn), &
               s_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
               sfields(id%salt)%off, sfields(id%salt)%ndims)
       end if
       s_bkginc(:,:,:) = s_bkginc(:,:,:) * tmask(:,:,:)

       if (id%temp>0 .or. id%salt>0) then
          ! Fill halo regions
          call lbc_lnk('distribute_state_pdaf', t_bkginc, 'T', &
               1., s_bkginc, 'T', 1.)
       end if

       ! U
       if (id%uvel > 0 .and. update_vel) then
          call state2field_inc(state_p, &
               uu(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, Nnn), &
               u_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
               sfields(id%uvel)%off, sfields(id%uvel)%ndims)
          u_bkginc(:,:,:) = u_bkginc(:,:,:) * umask(:,:,:)
       end if

       ! V
       if (id%vvel > 0 .and. update_vel) then
          call state2field_inc(state_p, &
               vv(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, Nnn), &
               v_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
               sfields(id%vvel)%off, sfields(id%vvel)%ndims)
          v_bkginc(:,:,:) = v_bkginc(:,:,:) * vmask(:,:,:)
       end if

       if (id%uvel>0 .or. id%vvel>0) then
          ! Fill halo regions
          call lbc_lnk('distribute_state_pdaf', u_bkginc, 'U', -1., v_bkginc, 'V', -1.)

          ! Apply divergence damping to the velocity increment
          call div_damping_filter(Nnn)
       end if


       ! For direct initialization we intialize the background arrays
       ! and apply the increments here. They are set the Euler flag
       if (ln_asmdin) then

           
          ! Ensure that the background arrays are allocated
          if (.not. allocated(ssh_bkg)) allocate(ssh_bkg(jpi,jpj))
          if (.not. allocated(t_bkg)) allocate(t_bkg(jpi,jpj,jpk))
          if (.not. allocated(s_bkg)) allocate(s_bkg(jpi,jpj,jpk))
          if (.not. allocated(u_bkg)) allocate(u_bkg(jpi,jpj,jpk))
          if (.not. allocated(v_bkg)) allocate(v_bkg(jpi,jpj,jpk))

          ! Initialize background fields for use in asm_inc routines
          ssh_bkg = ssh(:,:, Nnn)
          t_bkg = ts(:,:,:,jp_tem, Nnn)
          s_bkg = ts(:,:,:,jp_sal, Nnn)
          u_bkg = uu(:,:,:,Nnn)
          v_bkg = vv(:,:,:,Nnn)
          ssh_bkg(:,:) = ssh_bkg(:,:) * tmask(:,:,1)
          t_bkg(:,:,:) = t_bkg(:,:,:) * tmask(:,:,:)
          s_bkg(:,:,:) = s_bkg(:,:,:) * tmask(:,:,:)
          u_bkg(:,:,:) = u_bkg(:,:,:) * umask(:,:,:)
          v_bkg(:,:,:) = v_bkg(:,:,:) * vmask(:,:,:)

          if (mype_ens==0) &
               write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Apply full increment in ASMINC'

          ! Apply assimilation increment
          IF( ln_trainc )   CALL tra_asm_inc( next_inc, Nbb, Nnn, ts    , Nrhs )      ! Tracers
          IF( ln_dyninc )   CALL dyn_asm_inc( next_inc, Nbb, Nnn, uu, vv, Nrhs )      ! Dynamics
          IF( ln_sshinc )   CALL ssh_asm_inc( next_inc, Nbb, Nnn )                    ! SSH

       else
          if (mype_ens==0) &
               write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Store increments for IAU'
       end if

    end if physics


! *************************************************
! *** Prepare increment array for BCG variables ***
! *************************************************

#if defined key_top_asm
    do i_bgcinc = 1, n_update_trc

       ! Get index in BGC tracer array
       i_tr = ids_update_trc(i_bgcinc)

       if (sv_trc(i_tr)) then

          ! Get field index in state vector
          id_var=id%trc(i_tr)

          if (sfields(id_var)%update) then

             ! Update prognostic BGC variables
             call state2field_inc(state_p, &
                  tr(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, i_tr, Nnn), &
                  bgc_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, i_bgcinc), &
                  sfields(id_var)%off, sfields(id_var)%ndims)

             ! Fill halo regions
             call lbc_lnk('distribute_state_pdaf', bgc_bkginc(:, :, :, i_bgcinc), 'T', &
                  1._pwp)

             bgc_bkginc(:,:,:,i_bgcinc) = bgc_bkginc(:,:,:,i_bgcinc) * tmask(:,:,:)
          end if
       end if
    end do
#endif

  end subroutine update_bkginc_pdaf


!-------------------------------------------------------------------------------
!< Routine to apply divergence damping
!!
!! The functionality inside the routine is copied from NEMO's
!! routine asm_inc_init
!!
  subroutine div_damping_filter(Kmm)

    use asminc, &
         only: u_bkginc, v_bkginc
    use par_oce, &
         only: jpkm1, Nis0, Nie0, Njs0, Nje0, nn_hls
#if defined key_qco || defined key_linssh
#else
    use dom_oce, &
         only: e1v, e2u, e3u, e3v, e3t, &
         r1_e1u, r1_e2v, umask, vmask
#endif
    use lbclnk, &
         only: lbc_lnk

    implicit none

! *** Arguments
    integer, intent(in) :: Kmm

! *** Local variables ***
    integer :: ji, jj, jt, jk                ! Counters
    real(pwp), allocatable ::   zhdiv(:,:)   ! 2D workspace
    integer :: i0, j0, ni_p, nj_p

! ***************************************
! *** Apply divergence damping filter ***
! ***************************************

#if defined key_qco || defined key_linssh
#else
    divdmp: if ( ln_dyninc .and. nn_divdmp > 0 ) then    ! Apply divergence damping filter

       ! Compute halo offset
       i0 = Nis0 - nn_hls
       j0 = Njs0 - nn_hls

       ! Local dimensions
       ni_p = Nie0 - Nis0 + 1
       nj_p = Nje0 - Njs0 + 1

       if (mype_ens==0) &
            write (*,'(a,4x,a,i4)') 'NEMO-PDAF', '--- apply divergence damping: nn_divdmp', nn_divdmp

       allocate( zhdiv(jpi,jpj) ) 

       do jt = 1, nn_divdmp

          do jk = 1, jpkm1           ! zhdiv = e1e1 * div
             zhdiv(:,:) = 0.0
             do jj = j0+1, nj_p+j0
                do ji = i0+1, ni_p+i0
                   zhdiv(ji,jj) = (  e2u(ji  ,jj) * e3u(ji  ,jj,jk,Kmm) * u_bkginc(ji  ,jj,jk)    &
                        &            - e2u(ji-1,jj) * e3u(ji-1,jj,jk,Kmm) * u_bkginc(ji-1,jj,jk)    &
                        &            + e1v(ji,jj  ) * e3v(ji,jj  ,jk,Kmm) * v_bkginc(ji,jj  ,jk)    &
                        &            - e1v(ji,jj-1) * e3v(ji,jj-1,jk,Kmm) * v_bkginc(ji,jj-1,jk)  ) &
                        &            / e3t(ji,jj,jk,Kmm)
                end do
             end do
             call lbc_lnk( 'asminc', zhdiv, 'T', 1. )   ! lateral boundary cond. (no sign change)

             do jj = j0+1, nj_p+j0
                do ji = i0+1, ni_p+i0
                   u_bkginc(ji,jj,jk) = u_bkginc(ji,jj,jk)                         &
                        + 0.2_pwp * ( zhdiv(ji+1,jj) - zhdiv(ji  ,jj) ) * r1_e1u(ji,jj) * umask(ji,jj,jk)
                   v_bkginc(ji,jj,jk) = v_bkginc(ji,jj,jk)                         &
                        + 0.2_pwp * ( zhdiv(ji,jj+1) - zhdiv(ji,jj  ) ) * r1_e2v(ji,jj) * vmask(ji,jj,jk) 
                end do
             end do
          end do

       end do

       deallocate( zhdiv ) 

    endif divdmp
#endif
  end subroutine div_damping_filter

!-------------------------------------------------------------------------------
!> Routine to deallocate the BGC increment array
!!
  subroutine init_iauweight(weights, nsteps, iaushape)

    implicit none

! *** Arguments ***
    real(pwp), intent(out) :: weights(:)  ! Array of IAU weights
    integer, intent(in) :: nsteps         ! number of steps
    integer, intent(in) :: iaushape       ! Shape of the IAU weight function

! *** Local variables
    integer :: i, imid       ! Counter
    real(pwp) :: wnorm       ! Weights norm
     
    weights(:) = 0._pwp

    if( shape_asmiau == 0 ) then

       wnorm = 1.0 / real(nsteps, pwp)
        
       ! Constant IAU forcing
       do i = 1, nsteps
          weights(i) = wnorm
       end do

    elseif (shape_asmiau == 1) then

       ! Linear hat-like, centred in middle of IAU interval 

       ! Compute the normalization factor
       wnorm = 0.0
       if( mod(nsteps, 2) == 0 ) then   ! Even number of time steps in IAU interval
          imid = nsteps / 2 
          do i = 1, imid
             wnorm = wnorm + real(i)
          end do
          wnorm = 2.0 * wnorm
       else                                ! Odd number of time steps in IAU interval
          imid = (nsteps + 1) / 2        
          do i = 1, imid - 1
             wnorm = wnorm + real(i)
          end do
          wnorm = 2.0 * wnorm + real( imid )
       endif
       wnorm = 1.0 / wnorm
       
       do i = 1, imid - 1
          weights(i) = real(i) * wnorm
       end do
       do i = imid, nsteps
          weights(i) = real(nsteps - i + 1) * wnorm
       end do
    endif

  end subroutine init_iauweight


!-------------------------------------------------------------------------------
!> Routine to deallocate the BGC increment array
!!
  subroutine asm_inc_deallocate_pdaf

    implicit none
     
#if defined key_top_asm
    if (allocated(bgc_bkginc)) deallocate(bgc_bkginc)
#endif

  end subroutine asm_inc_deallocate_pdaf

end module asminc_pdaf
