!> Distributing the state vector variables at initial time
!!
!! This routine initializes the full fields of the
!! model from the statevector of PDAF (first timestep).
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! **Calling Sequence**
!!
!!  - Called from: `PDAF_get_state` (as U_dist_state)
!!
subroutine distribute_state_init_pdaf(dim_p, state_p)

  use mod_kind_pdaf
  use parallel_pdaf, &
       only: mype_ens
  use statevector_pdaf, &
       only: sfields, id
  use nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, neuler, lbc_lnk, lbc_lnk_multi, &
       sshb, tsb, ub, vb, &
       sshn, tsn, un, vn
  use assimilation_pdaf, &
       only: ens_restart
  use transforms_pdaf, &
       only: state2field, transform_field_mv
#if defined key_top
  use statevector_pdaf, &
       only: sv_trc
  use nemo_pdaf, &
       only: jptra, trb, trn
#endif

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: i            ! Counter
  integer :: id_var       ! Index
  integer :: verbose      ! Control verbosity


! ***********************************
! *** Apply field transformations ***
! ***********************************

  if (mype_ens==0) then
     verbose = 1
  else
     verbose = 0
  end if

  call transform_field_mv(2, state_p, 0, verbose) !21


! ********************************************
! *** Distribute state vector 2d variables ***
! ********************************************

  ! Note: The loop limits account for the halo offsets i0 and j0

  coldstart: if (.not. ens_restart) then

     if (verbose==1) write (*,'(a,4x,a)') 'NEMO-PDAF', 'distribute state at initial time'

     ! SSH
     if (id%ssh > 0) then
        call state2field(state_p, sshn(1+i0:ni_p+i0, 1+j0:nj_p+j0), &
             sfields(id%ssh)%off, sfields(id%ssh)%ndims)

        ! Fill halo regions
        call lbc_lnk('distribute_state_pdaf', sshn, 'T', 1.)

        ! Update before field
        sshb = sshn
     endif


! ************************************
! Distribute state vector 3d variables
! ************************************

     ! T
     if (id%temp > 0) then
        call state2field(state_p, &
             tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
             sfields(id%temp)%off, sfields(id%temp)%ndims)
     end if

     ! S
     if (id%salt > 0) then
        call state2field(state_p, &
             tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_sal), &
             sfields(id%salt)%off, sfields(id%salt)%ndims)
     end if

     if (id%temp>0 .or. id%salt>0) then
        ! Fill halo regions
        call lbc_lnk_multi('distribute_state_pdaf', tsn(:, :, :, jp_tem), 'T', &
             1., tsn(:, :, :, jp_sal), 'T', 1.)

        ! Update before fields
        tsb(:,:,:,jp_tem) = tsn(:,:,:,jp_tem)
        tsb(:,:,:,jp_sal) = tsn(:,:,:,jp_sal)
     end if

     ! U
     if (id%uvel > 0) then
        call state2field(state_p, &
             un(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
             sfields(id%uvel)%off, sfields(id%uvel)%ndims)
     end if

     ! V
     if (id%vvel > 0) then
        call state2field(state_p, &
             vn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
             sfields(id%vvel)%off, sfields(id%vvel)%ndims)
     end if

     if (id%uvel>0 .or. id%vvel>0) then
        ! Fill halo regions
        call lbc_lnk_multi('distribute_state_pdaf', un, 'U', -1., vn, 'V', -1.)

        ! Update before fields
        ub = un
        vb = vn
     end if

#if defined key_top
     ! BGC
     do i = 1, jptra
        if (sv_trc(i)) then
           id_var=id%trc(i)
           call state2field(state_p, &
                trn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, sfields(id_var)%id_trc), &
                sfields(id_var)%off, sfields(id_var)%ndims)

           ! Fill halo regions
           call lbc_lnk('distribute_state_pdaf', trn(:, :, :, sfields(id_var)%id_trc), 'T', 1._pwp)
           trb(:, :, :, sfields(id_var)%id_trc) = trn(:, :, :, sfields(id_var)%id_trc)
        end if
     end do
#endif

     ! Set Euler step
     neuler = 0

  else coldstart

     ! Ensemble restart using restart fields from NEMO

     if (verbose==1) &
          write (*,'(a,4x,a)') 'NEMO-PDAF', 'Ensemble restart - distribute_state inactive'

  end if coldstart

end subroutine distribute_state_init_pdaf
