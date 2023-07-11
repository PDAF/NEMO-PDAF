!> Set dimension of local model state
!!
!! The routine is called by PDAF during the
!! analysis step in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model state on the current analysis
!! domain. In addition, the coordinates of this
!! domain are stored and the index arrays for the
!! local state vector and the mapping between global
!! to local state vectors are initialized.
!! 
!! - Called from: `PDAFomi_assimilate_local`/`assimilation_pdaf`
!!
subroutine init_dim_l_pdaf(step, domain_p, dim_l)

  use mod_kind_pdaf
  use assimilation_pdaf, &
       only: domain_coords, dim_state_p, id_lstate_in_pstate
  use statevector_pdaf, &
       only: n_fields, sfields, sfields_l
  use nemo_pdaf, &
       only: lons, lats, use_wet_state, nwet, wet_pts, &
       sdim2d, deg2rad

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step     !< Current time step
  integer, intent(in)  :: domain_p !< Current local analysis domain
  integer, intent(out) :: dim_l    !< Local state dimension

! *** Local variables
  integer :: i, cnt, ifield        ! Counters
  integer :: id_surf               ! state vector index of surface grid point
  integer :: id_i, id_j            ! Grid coordinates for local analysis domain


  ! ****************************************
  ! *** Initialize local state dimension ***
  ! ****************************************

  ! dimension = (number of 2D variables) 
  !     + (number of 3D variables * number of wet layers)

  if (allocated(sfields_l)) deallocate(sfields_l)
  allocate(sfields_l(n_fields))

  dim_l = 0
  do i = 1, n_fields

     sfields_l(i)%off = dim_l

     if (sfields(i)%ndims==2) then
        sfields_l(i)%dim = 1
     else
        sfields_l(i)%dim = wet_pts(3, domain_p)
     end if

     dim_l = dim_l + sfields_l(i)%dim
  end do


  ! **********************************************
  ! *** Initialize coordinates of local domain ***
  ! **********************************************

  ! get grid point indices
  id_i = wet_pts(6, domain_p)
  id_j = wet_pts(7, domain_p)

  ! Use T-values to get local coordinates
  ! the coordinates are stored in radians (as required by PDAFOMI)
  domain_coords(1) = lons(id_i, id_j) * deg2rad
  domain_coords(2) = lats(id_i, id_j) * deg2rad


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  if (allocated(id_lstate_in_pstate)) deallocate(id_lstate_in_pstate)
  allocate(id_lstate_in_pstate(dim_l))

  cnt = 0

  do ifield = 1, n_fields

     ! Set indices
     if (use_wet_state==1) then
        ! state vector contains full columns of surface wet points

        id_surf = domain_p + sfields(ifield)%off

        if (sfields(ifield)%ndims==3) then
           do i = 1, wet_pts(3,domain_p)
              cnt = cnt + 1
              id_lstate_in_pstate(cnt) = id_surf + (i-1)*nwet
           enddo
        else
           cnt = cnt + 1
           id_lstate_in_pstate(cnt) = id_surf
        end if

     elseif (use_wet_state==2) then
        ! state vector only contains wet points - stored in leading vertical order

        if (sfields(ifield)%ndims==2) then
           cnt = cnt  + 1
           id_lstate_in_pstate(cnt) = domain_p + sfields(ifield)%off
        else
           id_surf = wet_pts(5, domain_p) + sfields(ifield)%off

           do i = 1, wet_pts(3,domain_p)
              cnt = cnt  + 1
              id_lstate_in_pstate(cnt) = id_surf + i - 1
           enddo
        end if

     else
        ! State vector contains full 3d box

        id_surf = wet_pts(4, domain_p) + sfields(ifield)%off

        if (sfields(ifield)%ndims==2) then
           cnt = cnt + 1
           id_lstate_in_pstate(cnt) = id_surf
        else
           do i = 1, wet_pts(3,domain_p)
              cnt = cnt + 1
              id_lstate_in_pstate(cnt) = id_surf + (i-1)*sdim2d
           enddo
        end if
     end if
  enddo

  if (id_lstate_in_pstate(wet_pts(3,domain_p)) > dim_state_p) then
     write(*,*) 'Error: please check the global indices for local state vector'
  endif

end subroutine init_dim_l_pdaf
