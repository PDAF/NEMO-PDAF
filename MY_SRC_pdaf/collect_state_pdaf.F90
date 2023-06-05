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
  use parallel_pdaf, &
       only: mype_model, task_id
  use statevector_pdaf, &
       only: sfields, id
  use nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, sshn, tsn, un, vn
#if defined key_top
  use statevector_pdaf, &
       only: n_trc, sv_trc
  use nemo_pdaf, &
       only: jptra, trn
#endif
  use transforms_pdaf, &
       only: field2state, transform_field_mv

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: i                ! Counter
  integer :: verbose          ! Control verbosity


  ! *********************************
  ! Collect state vector 2d variables
  ! *********************************

  ! Note: The calls account for the halo offsets i0 and j0

  ! SSH
  if (id%ssh > 0) then
     call field2state(sshn(1+i0:ni_p+i0, 1+j0:nj_p+j0), &
          state_p, &
          sfields(id%ssh)%off, sfields(id%ssh)%ndims)
  end if


  ! *********************************
  ! Collect state vector 3d variables
  ! *********************************

  ! T
  if (id%temp > 0) then
     call field2state(tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
          state_p, &
          sfields(id%temp)%off, sfields(id%temp)%ndims)
  end if

  ! S
  if (id%salt > 0) then
     call field2state(tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_sal), &
          state_p, &
          sfields(id%salt)%off, sfields(id%salt)%ndims)
  end if

  ! U
  if (id%uvel > 0) then
       call field2state(un(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
            state_p, &
            sfields(id%uvel)%off, sfields(id%uvel)%ndims)
  end if

  ! V
  if (id%vvel > 0) then
      call field2state(vn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
           state_p, &
           sfields(id%vvel)%off, sfields(id%vvel)%ndims)
  end if

#if defined key_top
  ! BGC
  do i = 1, jptra
     if (sv_trc(i)) then
        call field2state(trn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, sfields(id%trc(i))%id_trc), &
             state_p, &
             sfields(id%trc(i))%off, sfields(id%trc(i))%ndims)
     end if
  end do
#endif

  ! Aply field transformations
  if (mype_model==0 .and. task_id==1) then
     verbose = 1
  else
     verbose = 0
  end if

  call transform_field_mv(1, state_p, 0, verbose)

end subroutine collect_state_pdaf
