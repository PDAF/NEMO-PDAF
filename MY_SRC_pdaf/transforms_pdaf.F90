!> Module providing transformations between fields and state vector
!!
!! Routines to transform state vector
!! - mapping between state vector and model field
!! - state transformations (log, etc.)
!! - apply limits to values to fields in state vector
!!
module transforms_pdaf

  ! Include dimension information for model grid
  use mod_kind_pdaf
  use nemo_pdaf, &
       only: nlvls=>jpk, nj_p, ni_p, nwet, wet_pts, use_wet_state, i0, j0
  use statevector_pdaf, &
       only: id, sfields, n_fields
  use parallel_pdaf, &
       only: mype=>mype_model, task_id

  implicit none
  save

  interface field2state_missval
    module procedure field2state_missval_2d
    module procedure field2state_missval_3d
    module procedure field2state_missval_4d
  end interface 

  interface field2state
    module procedure field2state_2d
    module procedure field2state_3d
    module procedure field2state_4d
  end interface 
  
  interface state2field
    module procedure state2field_2d
    module procedure state2field_3d
    module procedure state2field_4d_dbl
    module procedure state2field_4d_sgl
    module procedure state2field_4d_dbl_mask
    module procedure state2field_4d_sgl_mask
  end interface state2field

  interface state2field_inc
     module procedure state2field_inc_2d
     module procedure state2field_inc_3d
  end interface

contains
!===============================================================================


!> Convert from NEMO model field to state vector
!!
  subroutine field2state_4d(field, state, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: field(:,:,:,:)   !< Model field
    real(pwp), intent(inout) :: state(:)         !< State vector
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize state vector from model field ***

    if (use_wet_state==1) then
       
       do k = 1, n_levels
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             if (k <= wet_pts(3,i)) then
                state(cnt) = field(wet_pts(6, i), wet_pts(7, i), k, 1)
             end if
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!$OMP PARALLEL DO PRIVATE (i, k, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                state(cnt + k) = field(wet_pts(6, i), wet_pts(7, i), k, 1)
             end do
          end do
!$OMP END PARALLEL DO

       else

!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             state(cnt) = field(wet_pts(6, i), wet_pts(7, i), 1, 1)
          end do
!$OMP END PARALLEL DO

       end if

    else
       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1,nj_p
             do i = 1, ni_p
                state(cnt) = field(i,j,k,1)
                cnt = cnt + 1
             enddo
          enddo
       enddo
    end if

  end subroutine field2state_4d
!==============================================================================

!> Convert from NEMO model field to state vector
!!
  subroutine field2state_3d(field, state, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: field(:,:,:)     !< Model field
    real(pwp), intent(inout) :: state(:)         !< State vector
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize state vector from model field ***

    if (use_wet_state==1) then
       
       do k = 1, n_levels
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             if (k <= wet_pts(3,i)) then
                state(cnt) = field(wet_pts(6, i), wet_pts(7, i), k)
             end if
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                state(cnt + k) = field(wet_pts(6, i), wet_pts(7, i), k)
             end do
          end do

       else
          do i = 1, nwet
             cnt = offset + i
             state(cnt) = field(wet_pts(6, i), wet_pts(7, i), 1)
          end do
       end if

    else
       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1,nj_p
             do i = 1, ni_p
                state(cnt) = field(i,j,k)
                cnt = cnt + 1
             enddo
          enddo
       enddo
    end if

  end subroutine field2state_3d
!==============================================================================

  !> Convert from NEMO model field to state vector
!!
  subroutine field2state_2d(field, state, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: field(:,:)       !< Model field
    real(pwp), intent(inout) :: state(:)         !< State vector
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j
    integer :: cnt

! *** Initialize state vector from model field ***

    if ((use_wet_state==1) .or. (use_wet_state==2)) then
       do i = 1, nwet
          cnt = i + offset
          state(cnt) = field(wet_pts(6, i), wet_pts(7, i))
       end do
    else
       cnt = 1 + offset
       do j = 1,nj_p
          do i = 1, ni_p
             state(cnt) = field(i,j)
             cnt = cnt + 1
          enddo
       enddo
    end if

  end subroutine field2state_2d
!===============================================================================


!> Convert from NEMO model field to state vector
!!
  subroutine field2state_missval_4d(field, state, offset, ndims, missval)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: field(:,:,:,:)   !< Model field
    real(pwp), intent(inout) :: state(:)         !< State vector
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field
    real(pwp), intent(in)    :: missval          !< missing value

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize state vector from model field ***

    if (use_wet_state==1) then
       
       do k = 1, n_levels
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             if (k <= wet_pts(3,i)) then
                state(cnt) = field(wet_pts(6, i), wet_pts(7, i), k, 1)
             end if
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!$OMP PARALLEL DO PRIVATE (i, k, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                state(cnt + k) = field(wet_pts(6, i), wet_pts(7, i), k, 1)
             end do
          end do
!$OMP END PARALLEL DO

       else

!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             state(cnt) = field(wet_pts(6, i), wet_pts(7, i), 1, 1)
          end do
!$OMP END PARALLEL DO

       end if

    else
       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1,nj_p
             do i = 1, ni_p
                if (abs(field(i,j,k,1)- missval) > 1e-10) then 
                   state(cnt) = field(i,j,k,1)
                else
                   state(cnt) = 0.0_pwp
                endif
                cnt = cnt + 1
             enddo
          enddo
       enddo
    end if

  end subroutine field2state_missval_4d
!==============================================================================

!> Convert from NEMO model field to state vector
!!
  subroutine field2state_missval_3d(field, state, offset, ndims, missval)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: field(:,:,:)     !< Model field
    real(pwp), intent(inout) :: state(:)         !< State vector
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field
    real(pwp), intent(in)    :: missval          !< missing value

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize state vector from model field ***

    if (use_wet_state==1) then
       
       do k = 1, n_levels
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             if (k <= wet_pts(3,i)) then
                state(cnt) = field(wet_pts(6, i), wet_pts(7, i), k)
             end if
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                state(cnt + k) = field(wet_pts(6, i), wet_pts(7, i), k)
             end do
          end do

       else
          do i = 1, nwet
             cnt = offset + i
             state(cnt) = field(wet_pts(6, i), wet_pts(7, i), 1)
          end do
       end if

    else
       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1,nj_p
             do i = 1, ni_p
                if (abs(field(i,j,k)- missval) > 1e-10) then 
                   state(cnt) = field(i,j,k)
                else
                   state(cnt) = 0.0_pwp
                endif
                cnt = cnt + 1
             enddo
          enddo
       enddo
    end if

  end subroutine field2state_missval_3d
!==============================================================================

  !> Convert from NEMO model field to state vector
!!
  subroutine field2state_missval_2d(field, state, offset, ndims, missval)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: field(:,:)       !< Model field
    real(pwp), intent(inout) :: state(:)         !< State vector
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field
    real(pwp), intent(in)    :: missval          !< missing value

! *** Local variables ***
    integer :: i, j
    integer :: cnt

! *** Initialize state vector from model field ***

    if ((use_wet_state==1) .or. (use_wet_state==2)) then
       do i = 1, nwet
          cnt = i + offset
          state(cnt) = field(wet_pts(6, i), wet_pts(7, i))
       end do
    else
       cnt = 1 + offset
       do j = 1,nj_p
          do i = 1, ni_p
             if (abs(field(i,j)- missval) > 1e-10) then 
                state(cnt) = field(i,j)
             else
                state(cnt) = 0.0_pwp
             endif
             cnt = cnt + 1
          enddo
       enddo
    end if

  end subroutine field2state_missval_2d

!==============================================================================

!> Convert from state vector to NEMO model field
!!
  subroutine state2field_4d_dbl(state, field, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(pwp), intent(out)   :: field(:,:,:,:)   !< Model field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize model field from state vector

    if (use_wet_state==1) then
       
       do k = 1, nlvls
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt)
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt + k)
             end do
!!!!$OMP END PARALLEL DO
          end do

       else

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             field(wet_pts(6, i), wet_pts(7, i), 1, 1) = state(cnt)
!!!!$OMP END PARALLEL DO
          end do

       end if
    else

       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1, nj_p
             do i = 1, ni_p
!                if (tmask(i + i0, j + j0, k) == 1.0_pwp) then
                   field(i,j,k,1) = state(cnt) !convert to NEMO format (ntimec,nlvls,nlats,nlons)
!                end if
                cnt = cnt + 1
             enddo
          enddo
       enddo

    end if

  end subroutine state2field_4d_dbl
!==============================================================================

!> Convert from state vector to NEMO model field
!!
  subroutine state2field_4d_sgl(state, field, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(4), intent(out)     :: field(:,:,:,:)   !< Model field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize model field from state vector

    if (use_wet_state==1) then
       
       do k = 1, nlvls
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt)
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt + k)
             end do
!!!!$OMP END PARALLEL DO
          end do

       else

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             field(wet_pts(6, i), wet_pts(7, i), 1, 1) = state(cnt)
!!!!$OMP END PARALLEL DO
          end do

       end if
    else

       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1, nj_p
             do i = 1, ni_p
!                if (tmask(i + i0, j + j0, k) == 1.0_pwp) then
                   field(i,j,k,1) = state(cnt) !convert to NEMO format (ntimec,nlvls,nlats,nlons)
!                end if
                cnt = cnt + 1
             enddo
          enddo
       enddo

    end if

  end subroutine state2field_4d_sgl

!==============================================================================

!> Convert from state vector to NEMO model field
!!
  subroutine state2field_4d_dbl_mask(state, field, offset, ndims, mask)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(pwp), intent(out)   :: field(:,:,:,:)   !< Model field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field
    real(pwp), intent(in)    :: mask(:,:,:)      !< Land mask

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize model field from state vector

    if (use_wet_state==1) then
       
       do k = 1, nlvls
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt)
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt + k)
             end do
!!!!$OMP END PARALLEL DO
          end do

       else

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             field(wet_pts(6, i), wet_pts(7, i), 1, 1) = state(cnt)
!!!!$OMP END PARALLEL DO
          end do

       end if
    else

       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1, nj_p
             do i = 1, ni_p
                if (mask(i + i0, j + j0, k) == 1.0_pwp) then
                   field(i,j,k,1) = state(cnt) !convert to NEMO format (ntimec,nlvls,nlats,nlons)
                end if
                cnt = cnt + 1
             enddo
          enddo
       enddo

    end if

  end subroutine state2field_4d_dbl_mask
!==============================================================================

!> Convert from state vector to NEMO model field
!!
  subroutine state2field_4d_sgl_mask(state, field, offset, ndims, mask)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(4), intent(out)     :: field(:,:,:,:)   !< Model field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field
    real(pwp), intent(in)    :: mask(:,:,:)      !< Land mask

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize model field from state vector

    if (use_wet_state==1) then
       
       do k = 1, nlvls
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt)
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt + k)
             end do
!!!!$OMP END PARALLEL DO
          end do

       else

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             field(wet_pts(6, i), wet_pts(7, i), 1, 1) = state(cnt)
!!!!$OMP END PARALLEL DO
          end do

       end if
    else

       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1, nj_p
             do i = 1, ni_p
                if (mask(i + i0, j + j0, k) == 1.0_pwp) then
                   field(i,j,k,1) = state(cnt) !convert to NEMO format (ntimec,nlvls,nlats,nlons)
                end if
                cnt = cnt + 1
             enddo
          enddo
       enddo

    end if

  end subroutine state2field_4d_sgl_mask
!==============================================================================

!> Convert from state vector to NEMO model field
!!
  subroutine state2field_3d(state, field, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(pwp), intent(out)   :: field(:,:,:)     !< Model field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize model field from state vector

    if (use_wet_state==1) then
       
       do k = 1, nlvls
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             field(wet_pts(6, i), wet_pts(7, i), k) = state(cnt)
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                field(wet_pts(6, i), wet_pts(7, i), k) = state(cnt + k)
             end do
!!!!$OMP END PARALLEL DO
          end do

       else

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             field(wet_pts(6, i), wet_pts(7, i), 1) = state(cnt)
!!!!$OMP END PARALLEL DO
          end do

       end if
    else

       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1, nj_p
             do i = 1, ni_p
                field(i,j,k) = state(cnt) !convert to NEMO format (ntimec,nlvls,nlats,nlons)
                cnt = cnt + 1
             enddo
          enddo
       enddo

    end if

  end subroutine state2field_3d
!==============================================================================
!> Convert from state vector to NEMO model field
!!
  subroutine state2field_2d(state, field, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(pwp), intent(out)   :: field(:,:)       !< Model field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j
    integer :: cnt


! *** Initialize model field from state vector

    if ((use_wet_state==1) .or. (use_wet_state==2)) then
       do i = 1, nwet
          cnt = i + offset
          field(wet_pts(6, i), wet_pts(7, i)) = state(cnt)
       end do
    else
       cnt = 1 + offset
       do j = 1, nj_p
          do i = 1, ni_p
                field(i,j) = state(cnt) !convert to NEMO format (ntimec,nlvls,nlats,nlons)
             cnt = cnt + 1
          enddo
       enddo
    end if

  end subroutine state2field_2d
!==============================================================================

!> Convert from state vector to NEMO model increment field
!!
  subroutine state2field_inc_3d(state, field, inc, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(pwp), intent(in)    :: field(:,:,:)     !< Model field
    real(pwp), intent(inout) :: inc(:,:,:)       !< Increment field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize model field from state vector

    if (use_wet_state==1) then
       
       do k = 1, nlvls
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             inc(wet_pts(6, i), wet_pts(7, i), k) &
                  = state(cnt) - field(wet_pts(6, i), wet_pts(7, i), k)
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                inc(wet_pts(6, i), wet_pts(7, i), k) &
                     = state(cnt + k) - field(wet_pts(6, i), wet_pts(7, i), k)
             end do
!!!!$OMP END PARALLEL DO
          end do

       else

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             inc(wet_pts(6, i), wet_pts(7, i), 1) &
                  = state(cnt) - field(wet_pts(6, i), wet_pts(7, i), 1)
!!!!$OMP END PARALLEL DO
          end do

       end if
    else

       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1, nj_p
             do i = 1, ni_p
                inc(i,j,k) = state(cnt) - field(i,j,k)
                cnt = cnt + 1
             enddo
          enddo
       enddo

    end if

  end subroutine state2field_inc_3d
!==============================================================================
!> Convert from state vector to NEMO model field
!!
  subroutine state2field_inc_2d(state, field, inc, offset, ndims)

    implicit none

! *** Arguments ***
    real(pwp), intent(in)    :: state(:)         !< State vector
    real(pwp), intent(in)    :: field(:,:)       !< Model field
    real(pwp), intent(inout) :: inc(:,:)         !< Increment field
    integer,   intent(in)    :: offset           !< Offset in state vector
    integer,   intent(in)    :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j
    integer :: cnt


! *** Initialize model field from state vector

    if ((use_wet_state==1) .or. (use_wet_state==2)) then
       do i = 1, nwet
          cnt = i + offset
          inc(wet_pts(6, i), wet_pts(7, i)) = state(cnt) - field(wet_pts(6, i), wet_pts(7, i))
       end do
    else
       cnt = 1 + offset
       do j = 1, nj_p
          do i = 1, ni_p
             inc(i,j) = state(cnt) - field(i,j)
             cnt = cnt + 1
          enddo
       enddo
    end if

  end subroutine state2field_inc_2d
!==============================================================================
!> Transform field, e.g. to log and back
!!
  subroutine transform_field(type, trafo, shift, state, dim, off, var, verbose)

    implicit none

    integer,         intent(in)    :: type     !< Direction transformation
    integer,         intent(in)    :: trafo    !< Type of transformation
    real(pwp),       intent(in)    :: shift    !< constant for shifting value in transformation
    real(pwp),       intent(inout) :: state(:) !< State vector
    integer,         intent(in)    :: dim      !< dimension of field in state vector
    integer,         intent(in)    :: off      !< Offset of field in state vector
    character(len=*),intent(in)    :: var      !< Name of variable
    integer,         intent(in)    :: verbose  !< (1) to write screen output


    if (type==1) then
       ! Transformation from NEMO value to transformed value

       select case (trafo)
       case(0)
!          write(*,*) 'No Transformation of variable ', trim(var)
       case(1)
          if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
               'NEMO-PDAF', '--- apply log-10 transformation to ', trim(var)
          state(off+1 : off+dim) = log10(state(off+1 : off+dim) + shift)
       case(2)
          if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
               'NEMO-PDAF', '--- apply ln transformation to', trim(var)
          state(off+1 : off+dim) = log(state(off+1 : off+dim) + shift)
       case(3)
          write(*,*) 'no transformation- box cox still needs to be implemented'
       case DEFAULT
!          write(*,*) 'No Transformation of variable ', trim(var)
       end select

    elseif (type==2) then
       ! Transformation from transformed value to NEMO value

       select case (trafo)
       case(0)
!          write(*,*) 'No Transformation of bio limit'
       case(1)
          if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
               'NEMO-PDAF', '--- revert log-10 transformation to ', trim(var)
          state(off+1 : off+dim) = (10.D0**state(off+1 : off+dim))-shift
       case(2)
          if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
               'NEMO-PDAF', '--- revert ln transformation to', trim(var)
          state(off+1 : off+dim) = (exp(state(off+1 : off+dim)))-shift
       case(3)
          write(*,*) 'no transformation- box cox still needs to be implemented'
       case DEFAULT
!          write(*,*) 'No Transformation of bio variable'
       end select

    end if

  end subroutine transform_field

!==============================================================================

!> Transform all fields, e.g. to log and back
!!
  subroutine transform_field_mv(type, state, limits, verbose)

    implicit none

! *** Arguments ***
    integer,   intent(in)    :: type     !< Direction of transformation
    real(pwp), intent(inout) :: state(:) !< State vector
    integer,   intent(in)    :: limits   !< Whether to also apply limits
             !< (0) no limits applied, 
             !< (11) limits applied to original variable for type 1
             !< (12) limits applied to transformed variable for type 1
             !< (21) limits applied to original variable for type 2
             !< (22) limits applied to transformed variable for type 2
    integer,   intent(in)    :: verbose  !< (1) to write screen output

! *** Local variables ***
    integer           :: i,j        ! Counters
    integer           :: cnt, cnt2  ! Counters
    integer           :: trafo      ! Type of transformation
    real(pwp)         :: shift      ! constant for shifting value in transformation
    integer           :: dim        ! dimension of field in state vector
    integer           :: off        ! Offset of field in state vector
    character(len=10) :: var        ! Name of variable
    integer           :: dolimit    ! Whether to apply a min/max limit
    real(pwp)         :: max_limit  ! Maximum limiting value
    real(pwp)         :: min_limit  ! Minimum limiting value


    ! *** Apply limits before performing transformation ***
    if (type==1) then

          ! Apply limits to original variable
          if (limits==11) call var_limits_mv(state, verbose)
    elseif (type==2) then

          ! Apply limits to transformed variable
          if (limits==22) call var_limits_mv(state, verbose)
    end if


    ! *** Apply field transformations ***

    if (verbose>0) write(*,'(a, 4x, a)') &
         'NEMO-PDAF', 'Apply field transformations'

    do i = 1, n_fields

       ! Initialize values from sfields
       trafo = sfields(i)%transform
       shift = sfields(i)%trafo_shift
       dim = sfields(i)%dim
       off = sfields(i)%off
       var = sfields(i)%variable
       dolimit = sfields(i)%limit
       min_limit = sfields(i)%min_limit
       max_limit = sfields(i)%max_limit

       
       ! *** Transformations
       if (type==1) then
          ! Transformation from NEMO value to transformed value

          select case (trafo)
          case(0)
!             if (verbose>0) write(*,'(a, 1x, a, 1x, a)') &
!                  'No transformation of variable ', trim(var)
          case(1)
             if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
                  'NEMO-PDAF', '--- apply log-10 transformation to ', trim(var)
             state(off+1 : off+dim) = log10(state(off+1 : off+dim) + shift)
          case(2)
             if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
                  'NEMO-PDAF', '--- apply ln transformation to', trim(var)
             cnt=0
             cnt2=0

             do j = off+1, off+dim 
                if (state(j)>0.0_pwp) then 
                   state(j) = log(state(j) + shift)
                   cnt2 = cnt2+1
                else
                   state(j) = -30.0_pwp
                   cnt = cnt+1
                end if
             end do

             if (verbose>0 .and. cnt>0) write (*,'(a,8x,a,1x,a,2i9)') &
                  'NEMO-PDAF','--- number of <=0, >0:', trim(var), cnt, cnt2

          case DEFAULT
             write(*,*) 'ERROR: no valid variable transformation selected, type', trafo
          end select

       elseif (type==2) then

          ! Transformation from transformed value to NEMO value

          select case (trafo)
          case(0)
!             if (verbose>0) write(*,'(a, 1x, a, 1x, a)') &
!                  'No transformation of variable ', trim(var)
          case(1)
             if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
                  'NEMO-PDAF', '--- revert log-10 transformation of ', trim(var)
             state(off+1 : off+dim) = (10.D0**state(off+1 : off+dim))-shift
          case(2)
             if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
                  'NEMO-PDAF', '--- revert ln transformation of ', trim(var)
             cnt=0
             cnt2=0
             do j = off+1, off+dim 
                if (state(j)>-30.0_pwp) then 
                   cnt2 = cnt2+1
                   state(j) = exp(state(j))-shift 
                else
                   state(j) = 0.0_pwp
                   cnt = cnt + 1
                end if
             end do

             if (verbose>0 .and. cnt>0) write (*,'(a,8x,a,1x,a,2i9)') &
                  'NEMO-PDAF','--- number of <=0, >0:', trim(var), cnt, cnt2

          case DEFAULT
             write(*,*) 'ERROR: no valid variable transformation selected, type', trafo
          end select

       end if

    end do

    ! *** Apply limits after performing transformation ***
    if (type==1) then

          ! Forward transform: Apply limits to transformed variable
          if (limits==12) call var_limits_mv(state, verbose)
    elseif (type==2) then

          ! Backwar transform: Apply limits to original variable
          if (limits==21) call var_limits_mv(state, verbose)
    end if

  end subroutine transform_field_mv

!==============================================================================

!> Apply limits to variables
!!
  subroutine var_limits_mv(state, verbose)

    implicit none

! *** Arguments ***
    real(pwp), intent(inout) :: state(:) !< State vector
    integer,   intent(in)    :: verbose  !< (1) to write screen output

! *** Local variables ***
    integer           :: i, j, cnt, cnt1  ! Counters
    integer           :: dim        ! dimension of field in state vector
    integer           :: off        ! Offset of field in state vector
    character(len=10) :: var        ! Name of variable
    integer           :: dolimit    ! Whether to apply a min/max limit
    real(pwp)         :: max_limit  ! Maximum limiting value
    real(pwp)         :: min_limit  ! Minimum limiting value


    if (verbose>0) write(*,'(a, 4x, a)') &
         'NEMO-PDAF', 'Apply limits to model fields'

    do i = 1, n_fields

       ! Initialize values from sfields
       dim = sfields(i)%dim
       off = sfields(i)%off
       var = sfields(i)%variable
       dolimit = sfields(i)%limit
       min_limit = sfields(i)%min_limit
       max_limit = sfields(i)%max_limit


       ! *** Apply Limits

       cnt = 0
       cnt1 = 0
       if (dolimit == 1) then

          if (verbose>0) write(*,'(a, 4x, a, es12.3, 1x, a, 1x, a)') &
               'NEMO-PDAF', '--- apply min. limit',min_limit,'to ', trim(var)

          ! Apply minimum limit
          do j = off+1, off+dim
             if (state(j) < min_limit) then
                state(j) = min_limit
                cnt = cnt + 1
             end if
          end do

          if (cnt>0 .and. verbose>0) &
               write(*,'(a, 8x, a, i8)') 'NEMO-PDAF', '--- number of affected values', cnt

       elseif (dolimit == 2) then

          if (verbose>0) write(*,'(a, 4x, a, es12.3, 1x, a, 1x, a)') &
               'NEMO-PDAF', '--- apply max. limit',max_limit,'to ', trim(var)

          ! Apply maximum limit
          do j = off+1, off+dim
             if (state(j) > max_limit) then
                state(j) = max_limit
                cnt = cnt + 1
             end if
          end do

          if (cnt>0 .and. verbose>0) &
               write(*,'(a, 8x, a, i8)') 'NEMO-PDAF', '--- number of affected values', cnt

       elseif (dolimit == 3) then

          if (verbose>0) write(*,'(a, 4x, a, 2es12.3, 1x, a, 1x, a)') &
               'NEMO-PDAF', '--- apply min/max limits of',min_limit,max_limit, 'to ', trim(var)

          ! Apply minimum and maximum limits
          do j = off+1, off+dim
             if (state(j) < min_limit) then
                state(j) = min_limit
                cnt1 = cnt1 + 1
             elseif (state(j) > max_limit) then
                state(j) = max_limit
                cnt = cnt + 1
             end if
          end do

          if ((cnt>0 .or. cnt1>0) .and. verbose>0) &
               write(*,'(a, 8x, a, 2i8)') 'NEMO-PDAF', '--- number of affected values', cnt1, cnt
       else

          if (verbose>0) write(*,'(a, 4x, a, 1x, a)') &
               'NEMO-PDAF', '--- no limit applied to ', trim(var)
       end if

    end do

  end subroutine var_limits_mv

end module transforms_pdaf
