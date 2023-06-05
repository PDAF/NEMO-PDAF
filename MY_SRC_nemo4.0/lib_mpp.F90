MODULE lib_mpp
   !!======================================================================
   !!                       ***  MODULE  lib_mpp  ***
   !! Ocean numerics:  massively parallel processing library
   !!=====================================================================
   !! History :  OPA  !  1994  (M. Guyon, J. Escobar, M. Imbard)  Original code
   !!            7.0  !  1997  (A.M. Treguier)  SHMEM additions
   !!            8.0  !  1998  (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI
   !!                 !  1998  (J.M. Molines) Open boundary conditions
   !!   NEMO     1.0  !  2003  (J.M. Molines, G. Madec)  F90, free form
   !!                 !  2003  (J.M. Molines) add mpp_ini_north(_3d,_2d)
   !!             -   !  2004  (R. Bourdalle Badie)  isend option in mpi
   !!                 !  2004  (J.M. Molines) minloc, maxloc
   !!             -   !  2005  (G. Madec, S. Masson)  npolj=5,6 F-point & ice cases
   !!             -   !  2005  (R. Redler) Replacement of MPI_COMM_WORLD except for MPI_Abort
   !!             -   !  2005  (R. Benshila, G. Madec)  add extra halo case
   !!             -   !  2008  (R. Benshila) add mpp_ini_ice
   !!            3.2  !  2009  (R. Benshila) SHMEM suppression, north fold in lbc_nfd
   !!            3.2  !  2009  (O. Marti)    add mpp_ini_znl
   !!            4.0  !  2011  (G. Madec)  move ctl_ routines from in_out_manager
   !!            3.5  !  2012  (S.Mocavero, I. Epicoco) Add mpp_lnk_bdy_3d/2d routines to optimize the BDY comm.
   !!            3.5  !  2013  (C. Ethe, G. Madec)  message passing arrays as local variables 
   !!            3.5  !  2013  (S.Mocavero, I.Epicoco - CMCC) north fold optimizations
   !!            3.6  !  2015  (O. Tintó and M. Castrillo - BSC) Added '_multiple' case for 2D lbc and max
   !!            4.0  !  2017  (G. Madec) automatique allocation of array argument (use any 3rd dimension)
   !!             -   !  2017  (G. Madec) create generic.h90 files to generate all lbc and north fold routines
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ctl_stop      : update momentum and tracer Kz from a tke scheme
   !!   ctl_warn      : initialization, namelist read, and parameters control
   !!   ctl_opn       : Open file and check if required file is available.
   !!   ctl_nam       : Prints informations when an error occurs while reading a namelist
   !!   get_unit      : give the index of an unused logical unit
   !!----------------------------------------------------------------------
#if   defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!   lib_mpp_alloc : allocate mpp arrays
   !!   mynode        : indentify the processor unit
   !!   mpp_lnk       : interface (defined in lbclnk) for message passing of 2d or 3d arrays (mpp_lnk_2d, mpp_lnk_3d)
   !!   mpp_lnk_icb   : interface for message passing of 2d arrays with extra halo for icebergs (mpp_lnk_2d_icb)
   !!   mpprecv       :
   !!   mppsend       :
   !!   mppscatter    :
   !!   mppgather     :
   !!   mpp_min       : generic interface for mppmin_int , mppmin_a_int , mppmin_real, mppmin_a_real
   !!   mpp_max       : generic interface for mppmax_int , mppmax_a_int , mppmax_real, mppmax_a_real
   !!   mpp_sum       : generic interface for mppsum_int , mppsum_a_int , mppsum_real, mppsum_a_real
   !!   mpp_minloc    :
   !!   mpp_maxloc    :
   !!   mppsync       :
   !!   mppstop       :
   !!   mpp_ini_north : initialisation of north fold
   !!   mpp_lbc_north_icb : alternative to mpp_nfd for extra outer halo with icebergs
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE lbcnfd         ! north fold treatment
   USE in_out_manager ! I/O manager
#if defined key_USE_PDAF
   USE parallel_pdaf, ONLY: init_parallel_pdaf
#endif

   IMPLICIT NONE
   PRIVATE

   INTERFACE mpp_nfd
      MODULE PROCEDURE   mpp_nfd_2d    , mpp_nfd_3d    , mpp_nfd_4d
      MODULE PROCEDURE   mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
   END INTERFACE

   ! Interface associated to the mpp_lnk_... routines is defined in lbclnk
   PUBLIC   mpp_lnk_2d    , mpp_lnk_3d    , mpp_lnk_4d
   PUBLIC   mpp_lnk_2d_ptr, mpp_lnk_3d_ptr, mpp_lnk_4d_ptr
   !
!!gm  this should be useless
   PUBLIC   mpp_nfd_2d    , mpp_nfd_3d    , mpp_nfd_4d
   PUBLIC   mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
!!gm end
   !
   PUBLIC   ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam
   PUBLIC   mynode, mppstop, mppsync, mpp_comm_free
   PUBLIC   mpp_ini_north
   PUBLIC   mpp_lnk_2d_icb
   PUBLIC   mpp_lbc_north_icb
   PUBLIC   mpp_min, mpp_max, mpp_sum, mpp_minloc, mpp_maxloc
   PUBLIC   mpp_delay_max, mpp_delay_sum, mpp_delay_rcv
   PUBLIC   mppscatter, mppgather
   PUBLIC   mpp_ini_znl
   PUBLIC   mppsend, mpprecv                          ! needed by TAM and ICB routines
   PUBLIC   mpp_lnk_bdy_2d, mpp_lnk_bdy_3d, mpp_lnk_bdy_4d
   
   !! * Interfaces
   !! define generic interface for these routine as they are called sometimes
   !! with scalar arguments instead of array arguments, which causes problems
   !! for the compilation on AIX system as well as NEC and SGI. Ok on COMPACQ
   INTERFACE mpp_min
      MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
   END INTERFACE
   INTERFACE mpp_max
      MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
   END INTERFACE
   INTERFACE mpp_sum
      MODULE PROCEDURE mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real,   &
         &             mppsum_realdd, mppsum_a_realdd
   END INTERFACE
   INTERFACE mpp_minloc
      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
   END INTERFACE
   INTERFACE mpp_maxloc
      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
   END INTERFACE

   !! ========================= !!
   !!  MPI  variable definition !!
   !! ========================= !!
!$AGRIF_DO_NOT_TREAT
   INCLUDE 'mpif.h'
!$AGRIF_END_DO_NOT_TREAT

   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .TRUE.    !: mpp flag

   INTEGER, PARAMETER         ::   nprocmax = 2**10   ! maximun dimension (required to be a power of 2)

   INTEGER, PUBLIC ::   mppsize        ! number of process
   INTEGER, PUBLIC ::   mpprank        ! process number  [ 0 - size-1 ]
!$AGRIF_DO_NOT_TREAT
   INTEGER, PUBLIC ::   mpi_comm_oce   ! opa local communicator
!$AGRIF_END_DO_NOT_TREAT

   INTEGER :: MPI_SUMDD

   ! variables used for zonal integration
   INTEGER, PUBLIC ::   ncomm_znl       !: communicator made by the processors on the same zonal average
   LOGICAL, PUBLIC ::   l_znl_root      !: True on the 'left'most processor on the same row
   INTEGER         ::   ngrp_znl        !  group ID for the znl processors
   INTEGER         ::   ndim_rank_znl   !  number of processors on the same zonal average
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::   nrank_znl  ! dimension ndim_rank_znl, number of the procs into the same znl domain

   ! North fold condition in mpp_mpi with jpni > 1 (PUBLIC for TAM)
   INTEGER, PUBLIC ::   ngrp_world        !: group ID for the world processors
   INTEGER, PUBLIC ::   ngrp_opa          !: group ID for the opa processors
   INTEGER, PUBLIC ::   ngrp_north        !: group ID for the northern processors (to be fold)
   INTEGER, PUBLIC ::   ncomm_north       !: communicator made by the processors belonging to ngrp_north
   INTEGER, PUBLIC ::   ndim_rank_north   !: number of 'sea' processor in the northern line (can be /= jpni !)
   INTEGER, PUBLIC ::   njmppmax          !: value of njmpp for the processors of the northern line
   INTEGER, PUBLIC ::   north_root        !: number (in the comm_opa) of proc 0 in the northern comm
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE, SAVE ::   nrank_north   !: dimension ndim_rank_north

   ! Type of send : standard, buffered, immediate
   CHARACTER(len=1), PUBLIC ::   cn_mpi_send        !: type od mpi send/recieve (S=standard, B=bsend, I=isend)
   LOGICAL         , PUBLIC ::   l_isend = .FALSE.  !: isend use indicator (T if cn_mpi_send='I')
   INTEGER         , PUBLIC ::   nn_buffer          !: size of the buffer in case of mpi_bsend

   ! Communications summary report
   CHARACTER(len=128), DIMENSION(:), ALLOCATABLE ::   crname_lbc                   !: names of lbc_lnk calling routines
   CHARACTER(len=128), DIMENSION(:), ALLOCATABLE ::   crname_glb                   !: names of global comm calling routines
   CHARACTER(len=128), DIMENSION(:), ALLOCATABLE ::   crname_dlg                   !: names of delayed global comm calling routines
   INTEGER, PUBLIC                               ::   ncom_stp = 0                 !: copy of time step # istp
   INTEGER, PUBLIC                               ::   ncom_fsbc = 1                !: copy of sbc time step # nn_fsbc
   INTEGER, PUBLIC                               ::   ncom_dttrc = 1               !: copy of top time step # nn_dttrc
   INTEGER, PUBLIC                               ::   ncom_freq                    !: frequency of comm diagnostic
   INTEGER, PUBLIC , DIMENSION(:,:), ALLOCATABLE ::   ncomm_sequence               !: size of communicated arrays (halos)
   INTEGER, PARAMETER, PUBLIC                    ::   ncom_rec_max = 5000          !: max number of communication record
   INTEGER, PUBLIC                               ::   n_sequence_lbc = 0           !: # of communicated arraysvia lbc
   INTEGER, PUBLIC                               ::   n_sequence_glb = 0           !: # of global communications
   INTEGER, PUBLIC                               ::   n_sequence_dlg = 0           !: # of delayed global communications
   INTEGER, PUBLIC                               ::   numcom = -1                  !: logical unit for communicaton report
   LOGICAL, PUBLIC                               ::   l_full_nf_update = .TRUE.    !: logical for a full (2lines) update of bc at North fold report
   INTEGER,                    PARAMETER, PUBLIC ::   nbdelay = 2       !: number of delayed operations
   !: name (used as id) of allreduce-delayed operations
   ! Warning: we must use the same character length in an array constructor (at least for gcc compiler)
   CHARACTER(len=32), DIMENSION(nbdelay), PUBLIC ::   c_delaylist = (/ 'cflice', 'fwb   ' /)
   !: component name where the allreduce-delayed operation is performed
   CHARACTER(len=3),  DIMENSION(nbdelay), PUBLIC ::   c_delaycpnt = (/ 'ICE'   , 'OCE' /)
   TYPE, PUBLIC ::   DELAYARR
      REAL(   wp), POINTER, DIMENSION(:) ::  z1d => NULL()
      COMPLEX(wp), POINTER, DIMENSION(:) ::  y1d => NULL()
   END TYPE DELAYARR
   TYPE( DELAYARR ), DIMENSION(nbdelay), PUBLIC, SAVE  ::   todelay         !: must have SAVE for default initialization of DELAYARR
   INTEGER,          DIMENSION(nbdelay), PUBLIC        ::   ndelayid = -1   !: mpi request id of the delayed operations

   ! timing summary report
   REAL(wp), DIMENSION(2), PUBLIC ::  waiting_time = 0._wp
   REAL(wp)              , PUBLIC ::  compute_time = 0._wp, elapsed_time = 0._wp
   
   REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::   tampon   ! buffer in case of bsend

   LOGICAL, PUBLIC ::   ln_nnogather                !: namelist control of northfold comms
   LOGICAL, PUBLIC ::   l_north_nogather = .FALSE.  !: internal control of northfold comms

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION mynode( ldtxt, ldname, kumnam_ref, kumnam_cfg, kumond, kstop, localComm )
      !!----------------------------------------------------------------------
      !!                  ***  routine mynode  ***
      !!
      !! ** Purpose :   Find processor unit
      !!----------------------------------------------------------------------
      CHARACTER(len=*),DIMENSION(:), INTENT(  out) ::   ldtxt        !
      CHARACTER(len=*)             , INTENT(in   ) ::   ldname       !
      INTEGER                      , INTENT(in   ) ::   kumnam_ref   ! logical unit for reference namelist
      INTEGER                      , INTENT(in   ) ::   kumnam_cfg   ! logical unit for configuration namelist
      INTEGER                      , INTENT(inout) ::   kumond       ! logical unit for namelist output
      INTEGER                      , INTENT(inout) ::   kstop        ! stop indicator
      INTEGER         , OPTIONAL   , INTENT(in   ) ::   localComm    !
      !
      INTEGER ::   mynode, ierr, code, ji, ii, ios
      LOGICAL ::   mpi_was_called
      !
      NAMELIST/nammpp/ cn_mpi_send, nn_buffer, jpni, jpnj, ln_nnogather
      !!----------------------------------------------------------------------
      !
      ii = 1
      WRITE(ldtxt(ii),*)                                                                  ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'mynode : mpi initialisation'                                    ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~ '                                                        ;   ii = ii + 1
      !
      REWIND( kumnam_ref )              ! Namelist nammpp in reference namelist: mpi variables
      READ  ( kumnam_ref, nammpp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nammpp in reference namelist', lwp )
      !
      REWIND( kumnam_cfg )              ! Namelist nammpp in configuration namelist: mpi variables
      READ  ( kumnam_cfg, nammpp, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nammpp in configuration namelist', lwp )
      !
      !                              ! control print
      WRITE(ldtxt(ii),*) '   Namelist nammpp'                                             ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      mpi send type          cn_mpi_send = ', cn_mpi_send       ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      size exported buffer   nn_buffer   = ', nn_buffer,' bytes';   ii = ii + 1
      !
      IF( jpni < 1 .OR. jpnj < 1  ) THEN
         WRITE(ldtxt(ii),*) '      jpni and jpnj will be calculated automatically' ;   ii = ii + 1
      ELSE
         WRITE(ldtxt(ii),*) '      processor grid extent in i         jpni = ',jpni       ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '      processor grid extent in j         jpnj = ',jpnj       ;   ii = ii + 1
      ENDIF

      WRITE(ldtxt(ii),*) '      avoid use of mpi_allgather at the north fold  ln_nnogather = ', ln_nnogather  ; ii = ii + 1

      CALL mpi_initialized ( mpi_was_called, code )
      IF( code /= MPI_SUCCESS ) THEN
         DO ji = 1, SIZE(ldtxt)
            IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
         END DO
         WRITE(*, cform_err)
         WRITE(*, *) 'lib_mpp: Error in routine mpi_initialized'
         CALL mpi_abort( mpi_comm_world, code, ierr )
      ENDIF

      IF( mpi_was_called ) THEN
         !
         SELECT CASE ( cn_mpi_send )
         CASE ( 'S' )                ! Standard mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Standard blocking mpi send (send)'             ;   ii = ii + 1
         CASE ( 'B' )                ! Buffer mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Buffer blocking mpi send (bsend)'              ;   ii = ii + 1
            IF( Agrif_Root() )   CALL mpi_init_oce( ldtxt, ii, ierr )
         CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
            WRITE(ldtxt(ii),*) '           Immediate non-blocking send (isend)'           ;   ii = ii + 1
            l_isend = .TRUE.
         CASE DEFAULT
            WRITE(ldtxt(ii),cform_err)                                                    ;   ii = ii + 1
            WRITE(ldtxt(ii),*) '           bad value for cn_mpi_send = ', cn_mpi_send     ;   ii = ii + 1
            kstop = kstop + 1
         END SELECT
         !
      ELSEIF ( PRESENT(localComm) .AND. .NOT. mpi_was_called ) THEN
         WRITE(ldtxt(ii),cform_err)                                                    ;   ii = ii + 1
         WRITE(ldtxt(ii),*) ' lib_mpp: You cannot provide a local communicator '          ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '          without calling MPI_Init before ! '                ;   ii = ii + 1
         kstop = kstop + 1
      ELSE
         SELECT CASE ( cn_mpi_send )
         CASE ( 'S' )                ! Standard mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Standard blocking mpi send (send)'             ;   ii = ii + 1
            CALL mpi_init( ierr )
         CASE ( 'B' )                ! Buffer mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Buffer blocking mpi send (bsend)'              ;   ii = ii + 1
            IF( Agrif_Root() )   CALL mpi_init_oce( ldtxt, ii, ierr )
         CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
            WRITE(ldtxt(ii),*) '           Immediate non-blocking send (isend)'           ;   ii = ii + 1
            l_isend = .TRUE.
            CALL mpi_init( ierr )
         CASE DEFAULT
            WRITE(ldtxt(ii),cform_err)                                                    ;   ii = ii + 1
            WRITE(ldtxt(ii),*) '           bad value for cn_mpi_send = ', cn_mpi_send     ;   ii = ii + 1
            kstop = kstop + 1
         END SELECT
         !
      ENDIF

      IF( PRESENT(localComm) ) THEN
         IF( Agrif_Root() ) THEN
            mpi_comm_oce = localComm
         ENDIF
      ELSE
         CALL mpi_comm_dup( mpi_comm_world, mpi_comm_oce, code)
         IF( code /= MPI_SUCCESS ) THEN
            DO ji = 1, SIZE(ldtxt)
               IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
            END DO
            WRITE(*, cform_err)
            WRITE(*, *) ' lib_mpp: Error in routine mpi_comm_dup'
            CALL mpi_abort( mpi_comm_world, code, ierr )
         ENDIF
      ENDIF

#if defined key_USE_PDAF
      ! Split the MPI communicator into MODEL, FILTER and COUPLE communicators for PDAF.
      CALL init_parallel_pdaf(mpi_comm_oce)
#endif

#if defined key_agrif
      IF( Agrif_Root() ) THEN
         CALL Agrif_MPI_Init(mpi_comm_oce)
      ELSE
         CALL Agrif_MPI_set_grid_comm(mpi_comm_oce)
      ENDIF
#endif

      CALL mpi_comm_rank( mpi_comm_oce, mpprank, ierr )
      CALL mpi_comm_size( mpi_comm_oce, mppsize, ierr )
      mynode = mpprank

      IF( mynode == 0 ) THEN
         CALL ctl_opn( kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. , 1 )
         WRITE(kumond, nammpp)      
      ENDIF
      !
      CALL MPI_OP_CREATE(DDPDD_MPI, .TRUE., MPI_SUMDD, ierr)
      !
   END FUNCTION mynode

   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_lnk_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_lnk_... routines
   !!                ptab   :   array or pointer of arrays on which the boundary condition is applied
   !!                cd_nat :   nature of array grid-points
   !!                psgn   :   sign used across the north fold boundary
   !!                kfld   :   optional, number of pt3d arrays
   !!                cd_mpp :   optional, fill the overlap area only
   !!                pval   :   optional, background value (used at closed boundaries)
   !!----------------------------------------------------------------------
   !
   !                       !==  2D array and array of 2D pointer  ==!
   !
#  define DIM_2d
#     define ROUTINE_LNK           mpp_lnk_2d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_2d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_LNK           mpp_lnk_3d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_3d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_LNK           mpp_lnk_4d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_4d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_4d

   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_nfd_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_nfd_... routines
   !!                ptab   :   array or pointer of arrays on which the boundary condition is applied
   !!                cd_nat :   nature of array grid-points
   !!                psgn   :   sign used across the north fold boundary
   !!                kfld   :   optional, number of pt3d arrays
   !!                cd_mpp :   optional, fill the overlap area only
   !!                pval   :   optional, background value (used at closed boundaries)
   !!----------------------------------------------------------------------
   !
   !                       !==  2D array and array of 2D pointer  ==!
   !
#  define DIM_2d
#     define ROUTINE_NFD           mpp_nfd_2d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_2d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_NFD           mpp_nfd_3d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_3d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_NFD           mpp_nfd_4d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_4d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_4d


   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_lnk_bdy_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_lnk_... routines
   !!                ptab   :   array or pointer of arrays on which the boundary condition is applied
   !!                cd_nat :   nature of array grid-points
   !!                psgn   :   sign used across the north fold boundary
   !!                kb_bdy :   BDY boundary set
   !!                kfld   :   optional, number of pt3d arrays
   !!----------------------------------------------------------------------
   !
   !                       !==  2D array and array of 2D pointer  ==!
   !
#  define DIM_2d
#     define ROUTINE_BDY           mpp_lnk_bdy_2d
#     include "mpp_bdy_generic.h90"
#     undef ROUTINE_BDY
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_BDY           mpp_lnk_bdy_3d
#     include "mpp_bdy_generic.h90"
#     undef ROUTINE_BDY
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_BDY           mpp_lnk_bdy_4d
#     include "mpp_bdy_generic.h90"
#     undef ROUTINE_BDY
#  undef DIM_4d

   !!----------------------------------------------------------------------
   !!
   !!   load_array  &   mpp_lnk_2d_9    à generaliser a 3D et 4D
   
   
   !!    mpp_lnk_sum_2d et 3D   ====>>>>>>   à virer du code !!!!
   
   
   !!----------------------------------------------------------------------



   SUBROUTINE mppsend( ktyp, pmess, kbytes, kdest, md_req )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppsend  ***
      !!
      !! ** Purpose :   Send messag passing array
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout) ::   pmess(*)   ! array of real
      INTEGER , INTENT(in   ) ::   kbytes     ! size of the array pmess
      INTEGER , INTENT(in   ) ::   kdest      ! receive process number
      INTEGER , INTENT(in   ) ::   ktyp       ! tag of the message
      INTEGER , INTENT(in   ) ::   md_req     ! argument for isend
      !!
      INTEGER ::   iflag
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( cn_mpi_send )
      CASE ( 'S' )                ! Standard mpi send (blocking)
         CALL mpi_send ( pmess, kbytes, mpi_double_precision, kdest , ktyp, mpi_comm_oce        , iflag )
      CASE ( 'B' )                ! Buffer mpi send (blocking)
         CALL mpi_bsend( pmess, kbytes, mpi_double_precision, kdest , ktyp, mpi_comm_oce        , iflag )
      CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
         ! be carefull, one more argument here : the mpi request identifier..
         CALL mpi_isend( pmess, kbytes, mpi_double_precision, kdest , ktyp, mpi_comm_oce, md_req, iflag )
      END SELECT
      !
   END SUBROUTINE mppsend


   SUBROUTINE mpprecv( ktyp, pmess, kbytes, ksource )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpprecv  ***
      !!
      !! ** Purpose :   Receive messag passing array
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout) ::   pmess(*)   ! array of real
      INTEGER , INTENT(in   ) ::   kbytes     ! suze of the array pmess
      INTEGER , INTENT(in   ) ::   ktyp       ! Tag of the recevied message
      INTEGER, OPTIONAL, INTENT(in) :: ksource    ! source process number
      !!
      INTEGER :: istatus(mpi_status_size)
      INTEGER :: iflag
      INTEGER :: use_source
      !!----------------------------------------------------------------------
      !
      ! If a specific process number has been passed to the receive call,
      ! use that one. Default is to use mpi_any_source
      use_source = mpi_any_source
      IF( PRESENT(ksource) )   use_source = ksource
      !
      CALL mpi_recv( pmess, kbytes, mpi_double_precision, use_source, ktyp, mpi_comm_oce, istatus, iflag )
      !
   END SUBROUTINE mpprecv


   SUBROUTINE mppgather( ptab, kp, pio )
      !!----------------------------------------------------------------------
      !!                   ***  routine mppgather  ***
      !!
      !! ** Purpose :   Transfert between a local subdomain array and a work
      !!     array which is distributed following the vertical level.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)      , INTENT(in   ) ::   ptab   ! subdomain input array
      INTEGER                           , INTENT(in   ) ::   kp     ! record length
      REAL(wp), DIMENSION(jpi,jpj,jpnij), INTENT(  out) ::   pio    ! subdomain input array
      !!
      INTEGER :: itaille, ierror   ! temporary integer
      !!---------------------------------------------------------------------
      !
      itaille = jpi * jpj
      CALL mpi_gather( ptab, itaille, mpi_double_precision, pio, itaille     ,   &
         &                            mpi_double_precision, kp , mpi_comm_oce, ierror )
      !
   END SUBROUTINE mppgather


   SUBROUTINE mppscatter( pio, kp, ptab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppscatter  ***
      !!
      !! ** Purpose :   Transfert between awork array which is distributed
      !!      following the vertical level and the local subdomain array.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpnij)  ::   pio    ! output array
      INTEGER                             ::   kp     ! Tag (not used with MPI
      REAL(wp), DIMENSION(jpi,jpj)        ::   ptab   ! subdomain array input
      !!
      INTEGER :: itaille, ierror   ! temporary integer
      !!---------------------------------------------------------------------
      !
      itaille = jpi * jpj
      !
      CALL mpi_scatter( pio, itaille, mpi_double_precision, ptab, itaille     ,   &
         &                            mpi_double_precision, kp  , mpi_comm_oce, ierror )
      !
   END SUBROUTINE mppscatter

   
   SUBROUTINE mpp_delay_sum( cdname, cdelay, y_in, pout, ldlast, kcom )
     !!----------------------------------------------------------------------
      !!                   ***  routine mpp_delay_sum  ***
      !!
      !! ** Purpose :   performed delayed mpp_sum, the result is received on next call
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   )               ::   cdname  ! name of the calling subroutine
      CHARACTER(len=*), INTENT(in   )               ::   cdelay  ! name (used as id) of the delayed operation
      COMPLEX(wp),      INTENT(in   ), DIMENSION(:) ::   y_in
      REAL(wp),         INTENT(  out), DIMENSION(:) ::   pout
      LOGICAL,          INTENT(in   )               ::   ldlast  ! true if this is the last time we call this routine
      INTEGER,          INTENT(in   ), OPTIONAL     ::   kcom
      !!
      INTEGER ::   ji, isz
      INTEGER ::   idvar
      INTEGER ::   ierr, ilocalcomm
      COMPLEX(wp), ALLOCATABLE, DIMENSION(:) ::   ytmp
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom

      isz = SIZE(y_in)
      
      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, ld_dlg = .TRUE. )

      idvar = -1
      DO ji = 1, nbdelay
         IF( TRIM(cdelay) == TRIM(c_delaylist(ji)) ) idvar = ji
      END DO
      IF ( idvar == -1 )   CALL ctl_stop( 'STOP',' mpp_delay_sum : please add a new delayed exchange for '//TRIM(cdname) )

      IF ( ndelayid(idvar) == 0 ) THEN         ! first call    with restart: %z1d defined in iom_delay_rst
         !                                       --------------------------
         IF ( SIZE(todelay(idvar)%z1d) /= isz ) THEN                  ! Check dimension coherence
            IF(lwp) WRITE(numout,*) ' WARNING: the nb of delayed variables in restart file is not the model one'
            DEALLOCATE(todelay(idvar)%z1d)
            ndelayid(idvar) = -1                                      ! do as if we had no restart
         ELSE
            ALLOCATE(todelay(idvar)%y1d(isz))
            todelay(idvar)%y1d(:) = CMPLX(todelay(idvar)%z1d(:), 0., wp)   ! create %y1d, complex variable needed by mpi_sumdd
         END IF
      ENDIF
      
      IF( ndelayid(idvar) == -1 ) THEN         ! first call without restart: define %y1d and %z1d from y_in with blocking allreduce
         !                                       --------------------------
         ALLOCATE(todelay(idvar)%z1d(isz), todelay(idvar)%y1d(isz))   ! allocate also %z1d as used for the restart
         CALL mpi_allreduce( y_in(:), todelay(idvar)%y1d(:), isz, MPI_DOUBLE_COMPLEX, mpi_sumdd, ilocalcomm, ierr )   ! get %y1d
         todelay(idvar)%z1d(:) = REAL(todelay(idvar)%y1d(:), wp)      ! define %z1d from %y1d
      ENDIF

      IF( ndelayid(idvar) > 0 )   CALL mpp_delay_rcv( idvar )         ! make sure %z1d is received

      ! send back pout from todelay(idvar)%z1d defined at previous call
      pout(:) = todelay(idvar)%z1d(:)

      ! send y_in into todelay(idvar)%y1d with a non-blocking communication
#if defined key_mpi2
      IF( ln_timing ) CALL tic_tac( .TRUE., ld_global = .TRUE.)
      CALL  mpi_allreduce( y_in(:), todelay(idvar)%y1d(:), isz, MPI_DOUBLE_COMPLEX, mpi_sumdd, ilocalcomm, ndelayid(idvar), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE., ld_global = .TRUE.)
#else
      CALL mpi_iallreduce( y_in(:), todelay(idvar)%y1d(:), isz, MPI_DOUBLE_COMPLEX, mpi_sumdd, ilocalcomm, ndelayid(idvar), ierr )
#endif

   END SUBROUTINE mpp_delay_sum

   
   SUBROUTINE mpp_delay_max( cdname, cdelay, p_in, pout, ldlast, kcom )
      !!----------------------------------------------------------------------
      !!                   ***  routine mpp_delay_max  ***
      !!
      !! ** Purpose :   performed delayed mpp_max, the result is received on next call
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   )                 ::   cdname  ! name of the calling subroutine
      CHARACTER(len=*), INTENT(in   )                 ::   cdelay  ! name (used as id) of the delayed operation
      REAL(wp),         INTENT(in   ), DIMENSION(:)   ::   p_in    ! 
      REAL(wp),         INTENT(  out), DIMENSION(:)   ::   pout    ! 
      LOGICAL,          INTENT(in   )                 ::   ldlast  ! true if this is the last time we call this routine
      INTEGER,          INTENT(in   ), OPTIONAL       ::   kcom
      !!
      INTEGER ::   ji, isz
      INTEGER ::   idvar
      INTEGER ::   ierr, ilocalcomm
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom

      isz = SIZE(p_in)

      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, ld_dlg = .TRUE. )

      idvar = -1
      DO ji = 1, nbdelay
         IF( TRIM(cdelay) == TRIM(c_delaylist(ji)) ) idvar = ji
      END DO
      IF ( idvar == -1 )   CALL ctl_stop( 'STOP',' mpp_delay_max : please add a new delayed exchange for '//TRIM(cdname) )

      IF ( ndelayid(idvar) == 0 ) THEN         ! first call    with restart: %z1d defined in iom_delay_rst
         !                                       --------------------------
         IF ( SIZE(todelay(idvar)%z1d) /= isz ) THEN                  ! Check dimension coherence
            IF(lwp) WRITE(numout,*) ' WARNING: the nb of delayed variables in restart file is not the model one'
            DEALLOCATE(todelay(idvar)%z1d)
            ndelayid(idvar) = -1                                      ! do as if we had no restart
         END IF
      ENDIF

      IF( ndelayid(idvar) == -1 ) THEN         ! first call without restart: define %z1d from p_in with a blocking allreduce
         !                                       --------------------------
         ALLOCATE(todelay(idvar)%z1d(isz))
         CALL mpi_allreduce( p_in(:), todelay(idvar)%z1d(:), isz, MPI_DOUBLE_PRECISION, mpi_max, ilocalcomm, ierr )   ! get %z1d
      ENDIF

      IF( ndelayid(idvar) > 0 )   CALL mpp_delay_rcv( idvar )         ! make sure %z1d is received

      ! send back pout from todelay(idvar)%z1d defined at previous call
      pout(:) = todelay(idvar)%z1d(:)

      ! send p_in into todelay(idvar)%z1d with a non-blocking communication
#if defined key_mpi2
      IF( ln_timing ) CALL tic_tac( .TRUE., ld_global = .TRUE.)
      CALL  mpi_allreduce( p_in(:), todelay(idvar)%z1d(:), isz, MPI_DOUBLE_PRECISION, mpi_max, ilocalcomm, ndelayid(idvar), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE., ld_global = .TRUE.)
#else
      CALL mpi_iallreduce( p_in(:), todelay(idvar)%z1d(:), isz, MPI_DOUBLE_PRECISION, mpi_max, ilocalcomm, ndelayid(idvar), ierr )
#endif

   END SUBROUTINE mpp_delay_max

   
   SUBROUTINE mpp_delay_rcv( kid )
      !!----------------------------------------------------------------------
      !!                   ***  routine mpp_delay_rcv  ***
      !!
      !! ** Purpose :  force barrier for delayed mpp (needed for restart) 
      !!
      !!----------------------------------------------------------------------
      INTEGER,INTENT(in   )      ::  kid 
      INTEGER ::   ierr
      !!----------------------------------------------------------------------
      IF( ndelayid(kid) /= -2 ) THEN  
#if ! defined key_mpi2
         IF( ln_timing ) CALL tic_tac( .TRUE., ld_global = .TRUE.)
         CALL mpi_wait( ndelayid(kid), MPI_STATUS_IGNORE, ierr )                        ! make sure todelay(kid) is received
         IF( ln_timing ) CALL tic_tac(.FALSE., ld_global = .TRUE.)
#endif
         IF( ASSOCIATED(todelay(kid)%y1d) )   todelay(kid)%z1d(:) = REAL(todelay(kid)%y1d(:), wp)  ! define %z1d from %y1d
         ndelayid(kid) = -2   ! add flag to know that mpi_wait was already called on kid
      ENDIF
   END SUBROUTINE mpp_delay_rcv

   
   !!----------------------------------------------------------------------
   !!    ***  mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real  ***
   !!   
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_MAX
#  define INTEGER_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmax_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmax_a_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef INTEGER_TYPE
!
#  define REAL_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmax_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmax_a_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef REAL_TYPE
#  undef OPERATION_MAX
   !!----------------------------------------------------------------------
   !!    ***  mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real  ***
   !!   
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_MIN
#  define INTEGER_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmin_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmin_a_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef INTEGER_TYPE
!
#  define REAL_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmin_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmin_a_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef REAL_TYPE
#  undef OPERATION_MIN

   !!----------------------------------------------------------------------
   !!    ***  mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real  ***
   !!   
   !!   Global sum of 1D array or a variable (integer, real or complex)
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_SUM
#  define INTEGER_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppsum_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppsum_a_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef INTEGER_TYPE
!
#  define REAL_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppsum_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppsum_a_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef REAL_TYPE
#  undef OPERATION_SUM

#  define OPERATION_SUM_DD
#  define COMPLEX_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppsum_realdd
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppsum_a_realdd
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef COMPLEX_TYPE
#  undef OPERATION_SUM_DD

   !!----------------------------------------------------------------------
   !!    ***  mpp_minloc2d, mpp_minloc3d, mpp_maxloc2d, mpp_maxloc3d
   !!   
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_MINLOC
#  define DIM_2d
#     define ROUTINE_LOC           mpp_minloc2d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_2d
#  define DIM_3d
#     define ROUTINE_LOC           mpp_minloc3d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_3d
#  undef OPERATION_MINLOC

#  define OPERATION_MAXLOC
#  define DIM_2d
#     define ROUTINE_LOC           mpp_maxloc2d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_2d
#  define DIM_3d
#     define ROUTINE_LOC           mpp_maxloc3d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_3d
#  undef OPERATION_MAXLOC

   SUBROUTINE mppsync()
      !!----------------------------------------------------------------------
      !!                  ***  routine mppsync  ***
      !!
      !! ** Purpose :   Massively parallel processors, synchroneous
      !!
      !!-----------------------------------------------------------------------
      INTEGER :: ierror
      !!-----------------------------------------------------------------------
      !
      CALL mpi_barrier( mpi_comm_oce, ierror )
      !
   END SUBROUTINE mppsync


   SUBROUTINE mppstop( ldfinal, ld_force_abort ) 
      !!----------------------------------------------------------------------
      !!                  ***  routine mppstop  ***
      !!
      !! ** purpose :   Stop massively parallel processors method
      !!
      !!----------------------------------------------------------------------
      LOGICAL, OPTIONAL, INTENT(in) :: ldfinal    ! source process number
      LOGICAL, OPTIONAL, INTENT(in) :: ld_force_abort    ! source process number
      LOGICAL ::   llfinal, ll_force_abort
      INTEGER ::   info
      !!----------------------------------------------------------------------
      llfinal = .FALSE.
      IF( PRESENT(ldfinal) ) llfinal = ldfinal
      ll_force_abort = .FALSE.
      IF( PRESENT(ld_force_abort) ) ll_force_abort = ld_force_abort
      !
      IF(ll_force_abort) THEN
         CALL mpi_abort( MPI_COMM_WORLD )
      ELSE
         CALL mppsync
         CALL mpi_finalize( info )
      ENDIF
      IF( .NOT. llfinal ) STOP 123
      !
   END SUBROUTINE mppstop


   SUBROUTINE mpp_comm_free( kcom )
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kcom
      !!
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      CALL MPI_COMM_FREE(kcom, ierr)
      !
   END SUBROUTINE mpp_comm_free


   SUBROUTINE mpp_ini_znl( kumout )
      !!----------------------------------------------------------------------
      !!               ***  routine mpp_ini_znl  ***
      !!
      !! ** Purpose :   Initialize special communicator for computing zonal sum
      !!
      !! ** Method  : - Look for processors in the same row
      !!              - Put their number in nrank_znl
      !!              - Create group for the znl processors
      !!              - Create a communicator for znl processors
      !!              - Determine if processor should write znl files
      !!
      !! ** output
      !!      ndim_rank_znl = number of processors on the same row
      !!      ngrp_znl = group ID for the znl processors
      !!      ncomm_znl = communicator for the ice procs.
      !!      n_znl_root = number (in the world) of proc 0 in the ice comm.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kumout   ! ocean.output logical units
      !
      INTEGER :: jproc      ! dummy loop integer
      INTEGER :: ierr, ii   ! local integer
      INTEGER, ALLOCATABLE, DIMENSION(:) ::   kwork
      !!----------------------------------------------------------------------
      !-$$     WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ngrp_world     : ', ngrp_world
      !-$$     WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - mpi_comm_world : ', mpi_comm_world
      !-$$     WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - mpi_comm_oce   : ', mpi_comm_oce
      !
      ALLOCATE( kwork(jpnij), STAT=ierr )
      IF( ierr /= 0 ) THEN
         WRITE(kumout, cform_err)
         WRITE(kumout,*) 'mpp_ini_znl : failed to allocate 1D array of length jpnij'
         CALL mppstop
      ENDIF

      IF( jpnj == 1 ) THEN
         ngrp_znl  = ngrp_world
         ncomm_znl = mpi_comm_oce
      ELSE
         !
         CALL MPI_ALLGATHER ( njmpp, 1, mpi_integer, kwork, 1, mpi_integer, mpi_comm_oce, ierr )
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - kwork pour njmpp : ', kwork
         !-$$        CALL flush(numout)
         !
         ! Count number of processors on the same row
         ndim_rank_znl = 0
         DO jproc=1,jpnij
            IF ( kwork(jproc) == njmpp ) THEN
               ndim_rank_znl = ndim_rank_znl + 1
            ENDIF
         END DO
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ndim_rank_znl : ', ndim_rank_znl
         !-$$        CALL flush(numout)
         ! Allocate the right size to nrank_znl
         IF (ALLOCATED (nrank_znl)) DEALLOCATE(nrank_znl)
         ALLOCATE(nrank_znl(ndim_rank_znl))
         ii = 0
         nrank_znl (:) = 0
         DO jproc=1,jpnij
            IF ( kwork(jproc) == njmpp) THEN
               ii = ii + 1
               nrank_znl(ii) = jproc -1
            ENDIF
         END DO
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - nrank_znl : ', nrank_znl
         !-$$        CALL flush(numout)

         ! Create the opa group
         CALL MPI_COMM_GROUP(mpi_comm_oce,ngrp_opa,ierr)
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ngrp_opa : ', ngrp_opa
         !-$$        CALL flush(numout)

         ! Create the znl group from the opa group
         CALL MPI_GROUP_INCL  ( ngrp_opa, ndim_rank_znl, nrank_znl, ngrp_znl, ierr )
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ngrp_znl ', ngrp_znl
         !-$$        CALL flush(numout)

         ! Create the znl communicator from the opa communicator, ie the pool of procs in the same row
         CALL MPI_COMM_CREATE ( mpi_comm_oce, ngrp_znl, ncomm_znl, ierr )
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ncomm_znl ', ncomm_znl
         !-$$        CALL flush(numout)
         !
      END IF

      ! Determines if processor if the first (starting from i=1) on the row
      IF ( jpni == 1 ) THEN
         l_znl_root = .TRUE.
      ELSE
         l_znl_root = .FALSE.
         kwork (1) = nimpp
         CALL mpp_min ( 'lib_mpp', kwork(1), kcom = ncomm_znl)
         IF ( nimpp == kwork(1)) l_znl_root = .TRUE.
      END IF

      DEALLOCATE(kwork)

   END SUBROUTINE mpp_ini_znl


   SUBROUTINE mpp_ini_north
      !!----------------------------------------------------------------------
      !!               ***  routine mpp_ini_north  ***
      !!
      !! ** Purpose :   Initialize special communicator for north folding
      !!      condition together with global variables needed in the mpp folding
      !!
      !! ** Method  : - Look for northern processors
      !!              - Put their number in nrank_north
      !!              - Create groups for the world processors and the north processors
      !!              - Create a communicator for northern processors
      !!
      !! ** output
      !!      njmppmax = njmpp for northern procs
      !!      ndim_rank_north = number of processors in the northern line
      !!      nrank_north (ndim_rank_north) = number  of the northern procs.
      !!      ngrp_world = group ID for the world processors
      !!      ngrp_north = group ID for the northern processors
      !!      ncomm_north = communicator for the northern procs.
      !!      north_root = number (in the world) of proc 0 in the northern comm.
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ierr
      INTEGER ::   jjproc
      INTEGER ::   ii, ji
      !!----------------------------------------------------------------------
      !
      njmppmax = MAXVAL( njmppt )
      !
      ! Look for how many procs on the northern boundary
      ndim_rank_north = 0
      DO jjproc = 1, jpnij
         IF( njmppt(jjproc) == njmppmax )   ndim_rank_north = ndim_rank_north + 1
      END DO
      !
      ! Allocate the right size to nrank_north
      IF (ALLOCATED (nrank_north)) DEALLOCATE(nrank_north)
      ALLOCATE( nrank_north(ndim_rank_north) )

      ! Fill the nrank_north array with proc. number of northern procs.
      ! Note : the rank start at 0 in MPI
      ii = 0
      DO ji = 1, jpnij
         IF ( njmppt(ji) == njmppmax   ) THEN
            ii=ii+1
            nrank_north(ii)=ji-1
         END IF
      END DO
      !
      ! create the world group
      CALL MPI_COMM_GROUP( mpi_comm_oce, ngrp_world, ierr )
      !
      ! Create the North group from the world group
      CALL MPI_GROUP_INCL( ngrp_world, ndim_rank_north, nrank_north, ngrp_north, ierr )
      !
      ! Create the North communicator , ie the pool of procs in the north group
      CALL MPI_COMM_CREATE( mpi_comm_oce, ngrp_north, ncomm_north, ierr )
      !
   END SUBROUTINE mpp_ini_north


   SUBROUTINE mpi_init_oce( ldtxt, ksft, code )
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_init.opa  ***
      !!
      !! ** Purpose :: export and attach a MPI buffer for bsend
      !!
      !! ** Method  :: define buffer size in namelist, if 0 no buffer attachment
      !!            but classical mpi_init
      !!
      !! History :: 01/11 :: IDRIS initial version for IBM only
      !!            08/04 :: R. Benshila, generalisation
      !!---------------------------------------------------------------------
      CHARACTER(len=*),DIMENSION(:), INTENT(  out) ::   ldtxt
      INTEGER                      , INTENT(inout) ::   ksft
      INTEGER                      , INTENT(  out) ::   code
      INTEGER                                      ::   ierr, ji
      LOGICAL                                      ::   mpi_was_called
      !!---------------------------------------------------------------------
      !
      CALL mpi_initialized( mpi_was_called, code )      ! MPI initialization
      IF ( code /= MPI_SUCCESS ) THEN
         DO ji = 1, SIZE(ldtxt)
            IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
         END DO
         WRITE(*, cform_err)
         WRITE(*, *) ' lib_mpp: Error in routine mpi_initialized'
         CALL mpi_abort( mpi_comm_world, code, ierr )
      ENDIF
      !
      IF( .NOT. mpi_was_called ) THEN
         CALL mpi_init( code )
         CALL mpi_comm_dup( mpi_comm_world, mpi_comm_oce, code )
         IF ( code /= MPI_SUCCESS ) THEN
            DO ji = 1, SIZE(ldtxt)
               IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
            END DO
            WRITE(*, cform_err)
            WRITE(*, *) ' lib_mpp: Error in routine mpi_comm_dup'
            CALL mpi_abort( mpi_comm_world, code, ierr )
         ENDIF
      ENDIF
      !
      IF( nn_buffer > 0 ) THEN
         WRITE(ldtxt(ksft),*) 'mpi_bsend, buffer allocation of  : ', nn_buffer   ;   ksft = ksft + 1
         ! Buffer allocation and attachment
         ALLOCATE( tampon(nn_buffer), stat = ierr )
         IF( ierr /= 0 ) THEN
            DO ji = 1, SIZE(ldtxt)
               IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
            END DO
            WRITE(*, cform_err)
            WRITE(*, *) ' lib_mpp: Error in ALLOCATE', ierr
            CALL mpi_abort( mpi_comm_world, code, ierr )
         END IF
         CALL mpi_buffer_attach( tampon, nn_buffer, code )
      ENDIF
      !
   END SUBROUTINE mpi_init_oce


   SUBROUTINE DDPDD_MPI( ydda, yddb, ilen, itype )
      !!---------------------------------------------------------------------
      !!   Routine DDPDD_MPI: used by reduction operator MPI_SUMDD
      !!
      !!   Modification of original codes written by David H. Bailey
      !!   This subroutine computes yddb(i) = ydda(i)+yddb(i)
      !!---------------------------------------------------------------------
      INTEGER                     , INTENT(in)    ::   ilen, itype
      COMPLEX(wp), DIMENSION(ilen), INTENT(in)    ::   ydda
      COMPLEX(wp), DIMENSION(ilen), INTENT(inout) ::   yddb
      !
      REAL(wp) :: zerr, zt1, zt2    ! local work variables
      INTEGER  :: ji, ztmp           ! local scalar
      !!---------------------------------------------------------------------
      !
      ztmp = itype   ! avoid compilation warning
      !
      DO ji=1,ilen
      ! Compute ydda + yddb using Knuth's trick.
         zt1  = real(ydda(ji)) + real(yddb(ji))
         zerr = zt1 - real(ydda(ji))
         zt2  = ((real(yddb(ji)) - zerr) + (real(ydda(ji)) - (zt1 - zerr))) &
                + aimag(ydda(ji)) + aimag(yddb(ji))

         ! The result is zt1 + zt2, after normalization.
         yddb(ji) = cmplx ( zt1 + zt2, zt2 - ((zt1 + zt2) - zt1),wp )
      END DO
      !
   END SUBROUTINE DDPDD_MPI


   SUBROUTINE mpp_lbc_north_icb( pt2d, cd_type, psgn, kextj)
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_icb  ***
      !!
      !! ** Purpose :   Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 and for 2d
      !!              array with outer extra halo
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4+kextj northern lines of the global domain on 1
      !!              processor and apply lbc north-fold on this sub array.
      !!              Then we scatter the north fold array back to the processors.
      !!              This routine accounts for an extra halo with icebergs
      !!              and assumes ghost rows and columns have been suppressed.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type  ! nature of pt3d grid-points
      !                                                     !   = T ,  U , V , F or W -points
      REAL(wp)                , INTENT(in   ) ::   psgn     ! = -1. the sign change across the
      !!                                                    ! north fold, =  1. otherwise
      INTEGER                 , INTENT(in   ) ::   kextj    ! Extra halo width at north fold
      !
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille, ildi, ilei, iilb
      INTEGER ::   ipj, ij, iproc
      !
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE  ::  ztab_e, znorthloc_e
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE  ::  znorthgloio_e
      !!----------------------------------------------------------------------
      !
      ipj=4
      ALLOCATE(        ztab_e(jpiglo, 1-kextj:ipj+kextj)       ,       &
     &            znorthloc_e(jpimax, 1-kextj:ipj+kextj)       ,       &
     &          znorthgloio_e(jpimax, 1-kextj:ipj+kextj,jpni)    )
      !
      ztab_e(:,:)      = 0._wp
      znorthloc_e(:,:) = 0._wp
      !
      ij = 1 - kextj
      ! put the last ipj+2*kextj lines of pt2d into znorthloc_e 
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         znorthloc_e(1:jpi,ij)=pt2d(1:jpi,jj)
         ij = ij + 1
      END DO
      !
      itaille = jpimax * ( ipj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER( znorthloc_e(1,1-kextj)    , itaille, MPI_DOUBLE_PRECISION,    &
         &                znorthgloio_e(1,1-kextj,1), itaille, MPI_DOUBLE_PRECISION,    &
         &                ncomm_north, ierr )
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      DO jr = 1, ndim_rank_north            ! recover the global north array
         iproc = nrank_north(jr) + 1
         ildi = nldit (iproc)
         ilei = nleit (iproc)
         iilb = nimppt(iproc)
         DO jj = 1-kextj, ipj+kextj
            DO ji = ildi, ilei
               ztab_e(ji+iilb-1,jj) = znorthgloio_e(ji,jj,jr)
            END DO
         END DO
      END DO

      ! 2. North-Fold boundary conditions
      ! ----------------------------------
      CALL lbc_nfd( ztab_e(:,1-kextj:ipj+kextj), cd_type, psgn, kextj )

      ij = 1 - kextj
      !! Scatter back to pt2d
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         DO ji= 1, jpi
            pt2d(ji,jj) = ztab_e(ji+nimpp-1,ij)
         END DO
         ij  = ij +1
      END DO
      !
      DEALLOCATE( ztab_e, znorthloc_e, znorthgloio_e )
      !
   END SUBROUTINE mpp_lbc_north_icb


   SUBROUTINE mpp_lnk_2d_icb( cdname, pt2d, cd_type, psgn, kexti, kextj )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_icb  ***
      !!
      !! ** Purpose :   Message passing management for 2d array (with extra halo for icebergs)
      !!                This routine receives a (1-kexti:jpi+kexti,1-kexti:jpj+kextj)
      !!                array (usually (0:jpi+1, 0:jpj+1)) from lbc_lnk_icb calls.
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    jpi    : first dimension of the local subdomain
      !!                    jpj    : second dimension of the local subdomain
      !!                    kexti  : number of columns for extra outer halo
      !!                    kextj  : number of rows for extra outer halo
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!----------------------------------------------------------------------
      CHARACTER(len=*)                                        , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      REAL(wp), DIMENSION(1-kexti:jpi+kexti,1-kextj:jpj+kextj), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)                                        , INTENT(in   ) ::   cd_type  ! nature of ptab array grid-points
      REAL(wp)                                                , INTENT(in   ) ::   psgn     ! sign used across the north fold
      INTEGER                                                 , INTENT(in   ) ::   kexti    ! extra i-halo width
      INTEGER                                                 , INTENT(in   ) ::   kextj    ! extra j-halo width
      !
      INTEGER  ::   jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! local integers
      INTEGER  ::   ipreci, iprecj             !   -       -
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for key_mpi_isend
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for key_mpi_isend
      !!
      REAL(wp), DIMENSION(1-kexti:jpi+kexti,nn_hls+kextj,2) ::   r2dns, r2dsn
      REAL(wp), DIMENSION(1-kextj:jpj+kextj,nn_hls+kexti,2) ::   r2dwe, r2dew
      !!----------------------------------------------------------------------

      ipreci = nn_hls + kexti      ! take into account outer extra 2D overlap area
      iprecj = nn_hls + kextj

      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, 1, 1, 1, ld_lbc = .TRUE. )

      ! 1. standard boundary treatment
      ! ------------------------------
      ! Order matters Here !!!!
      !
      !                                      ! East-West boundaries
      !                                           !* Cyclic east-west
      IF( l_Iperio ) THEN
         pt2d(1-kexti:     1   ,:) = pt2d(jpim1-kexti: jpim1 ,:)       ! east
         pt2d(  jpi  :jpi+kexti,:) = pt2d(     2     :2+kexti,:)       ! west
         !
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(  1-kexti   :nn_hls   ,:) = 0._wp    ! east except at F-point
                                      pt2d(jpi-nn_hls+1:jpi+kexti,:) = 0._wp    ! west
      ENDIF
      !                                      ! North-South boundaries
      IF( l_Jperio ) THEN                         !* cyclic (only with no mpp j-split)
         pt2d(:,1-kextj:     1   ) = pt2d(:,jpjm1-kextj:  jpjm1)       ! north
         pt2d(:,  jpj  :jpj+kextj) = pt2d(:,     2     :2+kextj)       ! south
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(:,  1-kextj   :nn_hls   ) = 0._wp    ! north except at F-point
                                      pt2d(:,jpj-nn_hls+1:jpj+kextj) = 0._wp    ! south
      ENDIF
      !

      ! north fold treatment
      ! -----------------------
      IF( npolj /= 0 ) THEN
         !
         SELECT CASE ( jpni )
                   CASE ( 1 )     ;   CALL lbc_nfd          ( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
                   CASE DEFAULT   ;   CALL mpp_lbc_north_icb( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
         END SELECT
         !
      ENDIF

      ! 2. East and west directions exchange
      ! ------------------------------------
      ! we play with the neigbours AND the row number because of the periodicity
      !
      SELECT CASE ( nbondi )      ! Read Dirichlet lateral conditions
      CASE ( -1, 0, 1 )                ! all exept 2 (i.e. close case)
         iihom = jpi-nreci-kexti
         DO jl = 1, ipreci
            r2dew(:,jl,1) = pt2d(nn_hls+jl,:)
            r2dwe(:,jl,1) = pt2d(iihom +jl,:)
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, r2dwe(1-kextj,1,1), imigr, noea, ml_req1 )
         CALL mpprecv( 1, r2dew(1-kextj,1,2), imigr, noea )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 1, r2dew(1-kextj,1,1), imigr, nowe, ml_req1 )
         CALL mppsend( 2, r2dwe(1-kextj,1,1), imigr, noea, ml_req2 )
         CALL mpprecv( 1, r2dew(1-kextj,1,2), imigr, noea )
         CALL mpprecv( 2, r2dwe(1-kextj,1,2), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 1, r2dew(1-kextj,1,1), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, r2dwe(1-kextj,1,2), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      iihom = jpi - nn_hls
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, ipreci
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
         END DO
      END SELECT


      ! 3. North and south directions
      ! -----------------------------
      ! always closed : we play only with the neigbours
      !
      IF( nbondj /= 2 ) THEN      ! Read Dirichlet lateral conditions
         ijhom = jpj-nrecj-kextj
         DO jl = 1, iprecj
            r2dsn(:,jl,1) = pt2d(:,ijhom +jl)
            r2dns(:,jl,1) = pt2d(:,nn_hls+jl)
         END DO
      ENDIF
      !
      !                           ! Migrations
      imigr = iprecj * ( jpi + 2*kexti )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, r2dsn(1-kexti,1,1), imigr, nono, ml_req1 )
         CALL mpprecv( 3, r2dns(1-kexti,1,2), imigr, nono )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 3, r2dns(1-kexti,1,1), imigr, noso, ml_req1 )
         CALL mppsend( 4, r2dsn(1-kexti,1,1), imigr, nono, ml_req2 )
         CALL mpprecv( 3, r2dns(1-kexti,1,2), imigr, nono )
         CALL mpprecv( 4, r2dsn(1-kexti,1,2), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 3, r2dns(1-kexti,1,1), imigr, noso, ml_req1 )
         CALL mpprecv( 4, r2dsn(1-kexti,1,2), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      ijhom = jpj - nn_hls
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, iprecj
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
         END DO
      END SELECT
      !
   END SUBROUTINE mpp_lnk_2d_icb


   SUBROUTINE mpp_report( cdname, kpk, kpl, kpf, ld_lbc, ld_glb, ld_dlg )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_report  ***
      !!
      !! ** Purpose :   report use of mpp routines per time-setp
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*),           INTENT(in   ) ::   cdname      ! name of the calling subroutine
      INTEGER         , OPTIONAL, INTENT(in   ) ::   kpk, kpl, kpf
      LOGICAL         , OPTIONAL, INTENT(in   ) ::   ld_lbc, ld_glb, ld_dlg
      !!
      CHARACTER(len=128)                        ::   ccountname  ! name of a subroutine to count communications
      LOGICAL ::   ll_lbc, ll_glb, ll_dlg
      INTEGER ::    ji,  jj,  jk,  jh, jf, jcount   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      ll_lbc = .FALSE.
      IF( PRESENT(ld_lbc) ) ll_lbc = ld_lbc
      ll_glb = .FALSE.
      IF( PRESENT(ld_glb) ) ll_glb = ld_glb
      ll_dlg = .FALSE.
      IF( PRESENT(ld_dlg) ) ll_dlg = ld_dlg
      !
      ! find the smallest common frequency: default = frequency product, if multiple, choose the larger of the 2 frequency
      IF( ncom_dttrc /= 1 )   CALL ctl_stop( 'STOP', 'mpp_report, ncom_dttrc /= 1 not coded...' ) 
      ncom_freq = ncom_fsbc
      !
      IF ( ncom_stp == nit000+ncom_freq ) THEN   ! avoid to count extra communications in potential initializations at nit000
         IF( ll_lbc ) THEN
            IF( .NOT. ALLOCATED(ncomm_sequence) ) ALLOCATE( ncomm_sequence(ncom_rec_max,2) )
            IF( .NOT. ALLOCATED(    crname_lbc) ) ALLOCATE(     crname_lbc(ncom_rec_max  ) )
            n_sequence_lbc = n_sequence_lbc + 1
            IF( n_sequence_lbc > ncom_rec_max ) CALL ctl_stop( 'STOP', 'lib_mpp, increase ncom_rec_max' )   ! deadlock
            crname_lbc(n_sequence_lbc) = cdname     ! keep the name of the calling routine
            ncomm_sequence(n_sequence_lbc,1) = kpk*kpl   ! size of 3rd and 4th dimensions
            ncomm_sequence(n_sequence_lbc,2) = kpf       ! number of arrays to be treated (multi)
         ENDIF
         IF( ll_glb ) THEN
            IF( .NOT. ALLOCATED(crname_glb) ) ALLOCATE( crname_glb(ncom_rec_max) )
            n_sequence_glb = n_sequence_glb + 1
            IF( n_sequence_glb > ncom_rec_max ) CALL ctl_stop( 'STOP', 'lib_mpp, increase ncom_rec_max' )   ! deadlock
            crname_glb(n_sequence_glb) = cdname     ! keep the name of the calling routine
         ENDIF
         IF( ll_dlg ) THEN
            IF( .NOT. ALLOCATED(crname_dlg) ) ALLOCATE( crname_dlg(ncom_rec_max) )
            n_sequence_dlg = n_sequence_dlg + 1
            IF( n_sequence_dlg > ncom_rec_max ) CALL ctl_stop( 'STOP', 'lib_mpp, increase ncom_rec_max' )   ! deadlock
            crname_dlg(n_sequence_dlg) = cdname     ! keep the name of the calling routine
         ENDIF
      ELSE IF ( ncom_stp == nit000+2*ncom_freq ) THEN
         CALL ctl_opn( numcom, 'communication_report.txt', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., narea )
         WRITE(numcom,*) ' '
         WRITE(numcom,*) ' ------------------------------------------------------------'
         WRITE(numcom,*) ' Communication pattern report (second oce+sbc+top time step):'
         WRITE(numcom,*) ' ------------------------------------------------------------'
         WRITE(numcom,*) ' '
         WRITE(numcom,'(A,I4)') ' Exchanged halos : ', n_sequence_lbc
         jj = 0; jk = 0; jf = 0; jh = 0
         DO ji = 1, n_sequence_lbc
            IF ( ncomm_sequence(ji,1) .GT. 1 ) jk = jk + 1
            IF ( ncomm_sequence(ji,2) .GT. 1 ) jf = jf + 1
            IF ( ncomm_sequence(ji,1) .GT. 1 .AND. ncomm_sequence(ji,2) .GT. 1 ) jj = jj + 1
            jh = MAX (jh, ncomm_sequence(ji,1)*ncomm_sequence(ji,2))
         END DO
         WRITE(numcom,'(A,I3)') ' 3D Exchanged halos : ', jk
         WRITE(numcom,'(A,I3)') ' Multi arrays exchanged halos : ', jf
         WRITE(numcom,'(A,I3)') '   from which 3D : ', jj
         WRITE(numcom,'(A,I10)') ' Array max size : ', jh*jpi*jpj
         WRITE(numcom,*) ' '
         WRITE(numcom,*) ' lbc_lnk called'
         DO ji = 1, n_sequence_lbc - 1
            IF ( crname_lbc(ji) /= 'already counted' ) THEN
               ccountname = crname_lbc(ji)
               crname_lbc(ji) = 'already counted'
               jcount = 1
               DO jj = ji + 1, n_sequence_lbc
                  IF ( ccountname ==  crname_lbc(jj) ) THEN
                     jcount = jcount + 1
                     crname_lbc(jj) = 'already counted'
                  END IF
               END DO
               WRITE(numcom,'(A, I4, A, A)') ' - ', jcount,' times by subroutine ', TRIM(ccountname)
            END IF
         END DO
         IF ( crname_lbc(n_sequence_lbc) /= 'already counted' ) THEN
            WRITE(numcom,'(A, I4, A, A)') ' - ', 1,' times by subroutine ', TRIM(crname_lbc(ncom_rec_max))
         END IF
         WRITE(numcom,*) ' '
         IF ( n_sequence_glb > 0 ) THEN
            WRITE(numcom,'(A,I4)') ' Global communications : ', n_sequence_glb
            jj = 1
            DO ji = 2, n_sequence_glb
               IF( crname_glb(ji-1) /= crname_glb(ji) ) THEN
                  WRITE(numcom,'(A, I4, A, A)') ' - ', jj,' times by subroutine ', TRIM(crname_glb(ji-1))
                  jj = 0
               END IF
               jj = jj + 1 
            END DO
            WRITE(numcom,'(A, I4, A, A)') ' - ', jj,' times by subroutine ', TRIM(crname_glb(n_sequence_glb))
            DEALLOCATE(crname_glb)
         ELSE
            WRITE(numcom,*) ' No MPI global communication '
         ENDIF
         WRITE(numcom,*) ' '
         IF ( n_sequence_dlg > 0 ) THEN
            WRITE(numcom,'(A,I4)') ' Delayed global communications : ', n_sequence_dlg
            jj = 1
            DO ji = 2, n_sequence_dlg
               IF( crname_dlg(ji-1) /= crname_dlg(ji) ) THEN
                  WRITE(numcom,'(A, I4, A, A)') ' - ', jj,' times by subroutine ', TRIM(crname_dlg(ji-1))
                  jj = 0
               END IF
               jj = jj + 1 
            END DO
            WRITE(numcom,'(A, I4, A, A)') ' - ', jj,' times by subroutine ', TRIM(crname_dlg(n_sequence_dlg))
            DEALLOCATE(crname_dlg)
         ELSE
            WRITE(numcom,*) ' No MPI delayed global communication '
         ENDIF
         WRITE(numcom,*) ' '
         WRITE(numcom,*) ' -----------------------------------------------'
         WRITE(numcom,*) ' '
         DEALLOCATE(ncomm_sequence)
         DEALLOCATE(crname_lbc)
      ENDIF
   END SUBROUTINE mpp_report

   
   SUBROUTINE tic_tac (ld_tic, ld_global)

    LOGICAL,           INTENT(IN) :: ld_tic
    LOGICAL, OPTIONAL, INTENT(IN) :: ld_global
    REAL(wp), DIMENSION(2), SAVE :: tic_wt
    REAL(wp),               SAVE :: tic_ct = 0._wp
    INTEGER :: ii

    IF( ncom_stp <= nit000 ) RETURN
    IF( ncom_stp == nitend ) RETURN
    ii = 1
    IF( PRESENT( ld_global ) ) THEN
       IF( ld_global ) ii = 2
    END IF
    
    IF ( ld_tic ) THEN
       tic_wt(ii) = MPI_Wtime()                                                    ! start count tic->tac (waiting time)
       IF ( tic_ct > 0.0_wp ) compute_time = compute_time + MPI_Wtime() - tic_ct   ! cumulate count tac->tic
    ELSE
       waiting_time(ii) = waiting_time(ii) + MPI_Wtime() - tic_wt(ii)              ! cumulate count tic->tac
       tic_ct = MPI_Wtime()                                                        ! start count tac->tic (waiting time)
    ENDIF
    
   END SUBROUTINE tic_tac

   
#else
   !!----------------------------------------------------------------------
   !!   Default case:            Dummy module        share memory computing
   !!----------------------------------------------------------------------
   USE in_out_manager

   INTERFACE mpp_sum
      MODULE PROCEDURE mppsum_int, mppsum_a_int, mppsum_real, mppsum_a_real, mppsum_realdd, mppsum_a_realdd
   END INTERFACE
   INTERFACE mpp_max
      MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
   END INTERFACE
   INTERFACE mpp_min
      MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
   END INTERFACE
   INTERFACE mpp_minloc
      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
   END INTERFACE
   INTERFACE mpp_maxloc
      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
   END INTERFACE

   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .FALSE.      !: mpp flag
   LOGICAL, PUBLIC            ::   ln_nnogather          !: namelist control of northfold comms (needed here in case "key_mpp_mpi" is not used)
   INTEGER, PUBLIC            ::   mpi_comm_oce          ! opa local communicator

   INTEGER, PARAMETER, PUBLIC               ::   nbdelay = 0   ! make sure we don't enter loops: DO ji = 1, nbdelay
   CHARACTER(len=32), DIMENSION(1), PUBLIC  ::   c_delaylist = 'empty'
   CHARACTER(len=32), DIMENSION(1), PUBLIC  ::   c_delaycpnt = 'empty'
   LOGICAL, PUBLIC                          ::   l_full_nf_update = .TRUE.
   TYPE ::   DELAYARR
      REAL(   wp), POINTER, DIMENSION(:) ::  z1d => NULL()
      COMPLEX(wp), POINTER, DIMENSION(:) ::  y1d => NULL()
   END TYPE DELAYARR
   TYPE( DELAYARR ), DIMENSION(1), PUBLIC  ::   todelay              
   INTEGER,  PUBLIC, DIMENSION(1)           ::   ndelayid = -1
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION lib_mpp_alloc(kumout)          ! Dummy function
      INTEGER, INTENT(in) ::   kumout
      lib_mpp_alloc = 0
   END FUNCTION lib_mpp_alloc

   FUNCTION mynode( ldtxt, ldname, kumnam_ref, knumnam_cfg,  kumond , kstop, localComm ) RESULT (function_value)
      INTEGER, OPTIONAL            , INTENT(in   ) ::   localComm
      CHARACTER(len=*),DIMENSION(:) ::   ldtxt
      CHARACTER(len=*) ::   ldname
      INTEGER ::   kumnam_ref, knumnam_cfg , kumond , kstop
      IF( PRESENT( localComm ) ) mpi_comm_oce = localComm
      function_value = 0
      IF( .FALSE. )   ldtxt(:) = 'never done'
      CALL ctl_opn( kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. , 1 )
   END FUNCTION mynode

   SUBROUTINE mppsync                       ! Dummy routine
   END SUBROUTINE mppsync

   !!----------------------------------------------------------------------
   !!    ***  mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real  ***
   !!   
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_MAX
#  define INTEGER_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmax_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmax_a_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef INTEGER_TYPE
!
#  define REAL_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmax_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmax_a_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef REAL_TYPE
#  undef OPERATION_MAX
   !!----------------------------------------------------------------------
   !!    ***  mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real  ***
   !!   
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_MIN
#  define INTEGER_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmin_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmin_a_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef INTEGER_TYPE
!
#  define REAL_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppmin_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppmin_a_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef REAL_TYPE
#  undef OPERATION_MIN

   !!----------------------------------------------------------------------
   !!    ***  mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real  ***
   !!   
   !!   Global sum of 1D array or a variable (integer, real or complex)
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_SUM
#  define INTEGER_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppsum_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppsum_a_int
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef INTEGER_TYPE
!
#  define REAL_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppsum_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppsum_a_real
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef REAL_TYPE
#  undef OPERATION_SUM

#  define OPERATION_SUM_DD
#  define COMPLEX_TYPE
#  define DIM_0d
#     define ROUTINE_ALLREDUCE           mppsum_realdd
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_0d
#  define DIM_1d
#     define ROUTINE_ALLREDUCE           mppsum_a_realdd
#     include "mpp_allreduce_generic.h90"
#     undef ROUTINE_ALLREDUCE
#  undef DIM_1d
#  undef COMPLEX_TYPE
#  undef OPERATION_SUM_DD

   !!----------------------------------------------------------------------
   !!    ***  mpp_minloc2d, mpp_minloc3d, mpp_maxloc2d, mpp_maxloc3d
   !!   
   !!----------------------------------------------------------------------
   !!
#  define OPERATION_MINLOC
#  define DIM_2d
#     define ROUTINE_LOC           mpp_minloc2d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_2d
#  define DIM_3d
#     define ROUTINE_LOC           mpp_minloc3d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_3d
#  undef OPERATION_MINLOC

#  define OPERATION_MAXLOC
#  define DIM_2d
#     define ROUTINE_LOC           mpp_maxloc2d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_2d
#  define DIM_3d
#     define ROUTINE_LOC           mpp_maxloc3d
#     include "mpp_loc_generic.h90"
#     undef ROUTINE_LOC
#  undef DIM_3d
#  undef OPERATION_MAXLOC

   SUBROUTINE mpp_delay_sum( cdname, cdelay, y_in, pout, ldlast, kcom )
      CHARACTER(len=*), INTENT(in   )               ::   cdname  ! name of the calling subroutine
      CHARACTER(len=*), INTENT(in   )               ::   cdelay  ! name (used as id) of the delayed operation
      COMPLEX(wp),      INTENT(in   ), DIMENSION(:) ::   y_in
      REAL(wp),         INTENT(  out), DIMENSION(:) ::   pout
      LOGICAL,          INTENT(in   )               ::   ldlast  ! true if this is the last time we call this routine
      INTEGER,          INTENT(in   ), OPTIONAL     ::   kcom
      !
      pout(:) = REAL(y_in(:), wp)
   END SUBROUTINE mpp_delay_sum

   SUBROUTINE mpp_delay_max( cdname, cdelay, p_in, pout, ldlast, kcom )
      CHARACTER(len=*), INTENT(in   )               ::   cdname  ! name of the calling subroutine
      CHARACTER(len=*), INTENT(in   )               ::   cdelay  ! name (used as id) of the delayed operation
      REAL(wp),         INTENT(in   ), DIMENSION(:) ::   p_in
      REAL(wp),         INTENT(  out), DIMENSION(:) ::   pout
      LOGICAL,          INTENT(in   )               ::   ldlast  ! true if this is the last time we call this routine
      INTEGER,          INTENT(in   ), OPTIONAL     ::   kcom
      !
      pout(:) = p_in(:)
   END SUBROUTINE mpp_delay_max

   SUBROUTINE mpp_delay_rcv( kid )
      INTEGER,INTENT(in   )      ::  kid 
      WRITE(*,*) 'mpp_delay_rcv: You should not have seen this print! error?', kid
   END SUBROUTINE mpp_delay_rcv
   
   SUBROUTINE mppstop( ldfinal, ld_force_abort )
      LOGICAL, OPTIONAL, INTENT(in) :: ldfinal    ! source process number
      LOGICAL, OPTIONAL, INTENT(in) :: ld_force_abort    ! source process number
      STOP      ! non MPP case, just stop the run
   END SUBROUTINE mppstop

   SUBROUTINE mpp_ini_znl( knum )
      INTEGER :: knum
      WRITE(*,*) 'mpp_ini_znl: You should not have seen this print! error?', knum
   END SUBROUTINE mpp_ini_znl

   SUBROUTINE mpp_comm_free( kcom )
      INTEGER :: kcom
      WRITE(*,*) 'mpp_comm_free: You should not have seen this print! error?', kcom
   END SUBROUTINE mpp_comm_free
   
#endif

   !!----------------------------------------------------------------------
   !!   All cases:         ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam   routines
   !!----------------------------------------------------------------------

   SUBROUTINE ctl_stop( cd1, cd2, cd3, cd4, cd5 ,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_opa  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the error number (nstop) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nstop = nstop + 1

      ! force to open ocean.output file
      IF( numout == 6 ) CALL ctl_opn( numout, 'ocean.output', 'APPEND', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
       
      WRITE(numout,cform_err)
      IF( PRESENT(cd1 ) )   WRITE(numout,*) TRIM(cd1)
      IF( PRESENT(cd2 ) )   WRITE(numout,*) TRIM(cd2)
      IF( PRESENT(cd3 ) )   WRITE(numout,*) TRIM(cd3)
      IF( PRESENT(cd4 ) )   WRITE(numout,*) TRIM(cd4)
      IF( PRESENT(cd5 ) )   WRITE(numout,*) TRIM(cd5)
      IF( PRESENT(cd6 ) )   WRITE(numout,*) TRIM(cd6)
      IF( PRESENT(cd7 ) )   WRITE(numout,*) TRIM(cd7)
      IF( PRESENT(cd8 ) )   WRITE(numout,*) TRIM(cd8)
      IF( PRESENT(cd9 ) )   WRITE(numout,*) TRIM(cd9)
      IF( PRESENT(cd10) )   WRITE(numout,*) TRIM(cd10)

                               CALL FLUSH(numout    )
      IF( numstp     /= -1 )   CALL FLUSH(numstp    )
      IF( numrun     /= -1 )   CALL FLUSH(numrun    )
      IF( numevo_ice /= -1 )   CALL FLUSH(numevo_ice)
      !
      IF( cd1 == 'STOP' ) THEN
         WRITE(numout,*)  'huge E-R-R-O-R : immediate stop'
         CALL mppstop(ld_force_abort = .true.)
      ENDIF
      !
   END SUBROUTINE ctl_stop


   SUBROUTINE ctl_warn( cd1, cd2, cd3, cd4, cd5,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_warn  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the warning number (nwarn) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nwarn = nwarn + 1
      IF(lwp) THEN
         WRITE(numout,cform_war)
         IF( PRESENT(cd1 ) ) WRITE(numout,*) TRIM(cd1)
         IF( PRESENT(cd2 ) ) WRITE(numout,*) TRIM(cd2)
         IF( PRESENT(cd3 ) ) WRITE(numout,*) TRIM(cd3)
         IF( PRESENT(cd4 ) ) WRITE(numout,*) TRIM(cd4)
         IF( PRESENT(cd5 ) ) WRITE(numout,*) TRIM(cd5)
         IF( PRESENT(cd6 ) ) WRITE(numout,*) TRIM(cd6)
         IF( PRESENT(cd7 ) ) WRITE(numout,*) TRIM(cd7)
         IF( PRESENT(cd8 ) ) WRITE(numout,*) TRIM(cd8)
         IF( PRESENT(cd9 ) ) WRITE(numout,*) TRIM(cd9)
         IF( PRESENT(cd10) ) WRITE(numout,*) TRIM(cd10)
      ENDIF
      CALL FLUSH(numout)
      !
   END SUBROUTINE ctl_warn


   SUBROUTINE ctl_opn( knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_opn  ***
      !!
      !! ** Purpose :   Open file and check if required file is available.
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(  out) ::   knum      ! logical unit to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdfile    ! file name to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdstat    ! disposition specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdform    ! formatting specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdacce    ! access specifier
      INTEGER          , INTENT(in   ) ::   klengh    ! record length
      INTEGER          , INTENT(in   ) ::   kout      ! number of logical units for write
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      INTEGER, OPTIONAL, INTENT(in   ) ::   karea     ! proc number
      !
      CHARACTER(len=80) ::   clfile
      INTEGER           ::   iost
      !!----------------------------------------------------------------------
      !
      ! adapt filename
      ! ----------------
      clfile = TRIM(cdfile)
      IF( PRESENT( karea ) ) THEN
         IF( karea > 1 )   WRITE(clfile, "(a,'_',i4.4)") TRIM(clfile), karea-1
      ENDIF
#if defined key_agrif
      IF( .NOT. Agrif_Root() )   clfile = TRIM(Agrif_CFixed())//'_'//TRIM(clfile)
      knum=Agrif_Get_Unit()
#else
      knum=get_unit()
#endif
      IF( TRIM(cdfile) == '/dev/null' )   clfile = TRIM(cdfile)   ! force the use of /dev/null
      !
      iost=0
      IF( cdacce(1:6) == 'DIRECT' )  THEN         ! cdacce has always more than 6 characters
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat, RECL=klengh         , ERR=100, IOSTAT=iost )
      ELSE IF( TRIM(cdstat) == 'APPEND' )  THEN   ! cdstat can have less than 6 characters
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS='UNKNOWN', POSITION='APPEND', ERR=100, IOSTAT=iost )
      ELSE
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat                      , ERR=100, IOSTAT=iost )
      ENDIF
      IF( iost /= 0 .AND. TRIM(clfile) == '/dev/null' ) &   ! for windows
         &  OPEN(UNIT=knum,FILE='NUL', FORM=cdform, ACCESS=cdacce, STATUS=cdstat                      , ERR=100, IOSTAT=iost )   
      IF( iost == 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*) '     file   : ', TRIM(clfile),' open ok'
            WRITE(kout,*) '     unit   = ', knum
            WRITE(kout,*) '     status = ', cdstat
            WRITE(kout,*) '     form   = ', cdform
            WRITE(kout,*) '     access = ', cdacce
            WRITE(kout,*)
         ENDIF
      ENDIF
100   CONTINUE
      IF( iost /= 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*)
            WRITE(kout,*) ' ===>>>> : bad opening file: ', TRIM(clfile)
            WRITE(kout,*) ' =======   ===  '
            WRITE(kout,*) '           unit   = ', knum
            WRITE(kout,*) '           status = ', cdstat
            WRITE(kout,*) '           form   = ', cdform
            WRITE(kout,*) '           access = ', cdacce
            WRITE(kout,*) '           iostat = ', iost
            WRITE(kout,*) '           we stop. verify the file '
            WRITE(kout,*)
         ELSE  !!! Force writing to make sure we get the information - at least once - in this violent STOP!!
            WRITE(*,*)
            WRITE(*,*) ' ===>>>> : bad opening file: ', TRIM(clfile)
            WRITE(*,*) ' =======   ===  '
            WRITE(*,*) '           unit   = ', knum
            WRITE(*,*) '           status = ', cdstat
            WRITE(*,*) '           form   = ', cdform
            WRITE(*,*) '           access = ', cdacce
            WRITE(*,*) '           iostat = ', iost
            WRITE(*,*) '           we stop. verify the file '
            WRITE(*,*)
         ENDIF
         CALL FLUSH( kout ) 
         STOP 'ctl_opn bad opening'
      ENDIF
      !
   END SUBROUTINE ctl_opn


   SUBROUTINE ctl_nam ( kios, cdnam, ldwp )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_nam  ***
      !!
      !! ** Purpose :   Informations when error while reading a namelist
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(inout) ::   kios    ! IO status after reading the namelist
      CHARACTER(len=*), INTENT(in   ) ::   cdnam   ! group name of namelist for which error occurs
      CHARACTER(len=5)                ::   clios   ! string to convert iostat in character for print
      LOGICAL         , INTENT(in   ) ::   ldwp    ! boolean term for print
      !!----------------------------------------------------------------------
      !
      WRITE (clios, '(I5.0)')   kios
      IF( kios < 0 ) THEN         
         CALL ctl_warn( 'end of record or file while reading namelist '   &
            &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      !
      IF( kios > 0 ) THEN
         CALL ctl_stop( 'misspelled variable in namelist '   &
            &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      kios = 0
      RETURN
      !
   END SUBROUTINE ctl_nam


   INTEGER FUNCTION get_unit()
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  get_unit  ***
      !!
      !! ** Purpose :   return the index of an unused logical unit
      !!----------------------------------------------------------------------
      LOGICAL :: llopn
      !!----------------------------------------------------------------------
      !
      get_unit = 15   ! choose a unit that is big enough then it is not already used in NEMO
      llopn = .TRUE.
      DO WHILE( (get_unit < 998) .AND. llopn )
         get_unit = get_unit + 1
         INQUIRE( unit = get_unit, opened = llopn )
      END DO
      IF( (get_unit == 999) .AND. llopn ) THEN
         CALL ctl_stop( 'get_unit: All logical units until 999 are used...' )
         get_unit = -1
      ENDIF
      !
   END FUNCTION get_unit

   !!----------------------------------------------------------------------
END MODULE lib_mpp
