subroutine get_nemo_indices(nemo_Nbb, nemo_Nnn, nemo_Nrhs)

#if defined key_qco   ||   defined key_linssh
  use stpmlf, &
       only: Nbb, Nnn, Nrhs
#else
  use step, &
       only: Nbb, Nnn, Nrhs
#endif

  implicit none

  integer, intent(out) :: nemo_Nbb
  integer, intent(out) :: nemo_Nnn
  integer, intent(out) :: nemo_Nrhs


  ! Hand over NEMO dimensions for use with nemo_pdaf
  
  nemo_Nbb = Nbb
  nemo_Nnn = Nnn
  nemo_Nrhs = Nrhs

end subroutine get_nemo_indices
