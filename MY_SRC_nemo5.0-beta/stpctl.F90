MODULE stpctl
   !!======================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!======================================================================
   !! History :  OPA  ! 1991-03  (G. Madec) Original code
   !!            6.0  ! 1992-06  (M. Imbard)
   !!            8.0  ! 1997-06  (A.M. Treguier)
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2009-07  (G. Madec)  Add statistic for time-spliting
   !!            3.7  ! 2016-09  (G. Madec)  Remove solver
   !!            4.0  ! 2017-04  (G. Madec)  regroup global communications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce ,  ONLY : ln_zad_Aimp       ! ocean vertical physics variables
   USE wet_dry,   ONLY : ll_wd, ssh_ref    ! reference depth for negative bathy
   !  
   USE diawri          ! Standard run outputs       (dia_wri_state routine)
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE eosbn2, ONLY: ln_SEOS, rn_b0
   USE timing          ! timing
   !
   USE netcdf          ! NetCDF library
#if defined key_USE_PDAF
   USE parallel_pdaf, ONLY : task_str ! Include ensemble task id
#endif
   USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_nan

   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_ctl           ! routine called by step.F90

   INTEGER, PARAMETER         ::   jpvar = 9
   INTEGER                    ::   nrunid   ! netcdf file id
   INTEGER, DIMENSION(jpvar)  ::   nvarid   ! netcdf variable id
   INTEGER, DIMENSION(3)      ::   j0oce    ! global integer wet point indicators
   LOGICAL                    ::   l0oce_T, l0oce_U, l0oce_V        ! Logical no wet point indicators
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: stpctl.F90 15023 2021-06-18 14:35:25Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_ctl( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Stop the run IF problem encountered by setting nstop > 0
      !!                Problems checked: |ssh| maximum larger than 10 m
      !!                                  |U|   maximum larger than 10 m/s 
      !!                                  negative sea surface salinity
      !!
      !! ** Actions :   "time.step" file = last ocean time-step
      !!                "run.stat"  file = run statistics
      !!                 nstop indicator sheared among all local domain
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER, INTENT(in   ) ::   Kmm      ! ocean time level index
      !!
      INTEGER, PARAMETER              ::   jptst = 5
      INTEGER                         ::   ji                                    ! dummy loop indices
      INTEGER                         ::   idtime, istatus
      INTEGER , DIMENSION(jptst)      ::   iareasum, iareamin, iareamax
      INTEGER , DIMENSION(3,jptst)    ::   iloc                                  ! min/max loc indices
      REAL(wp)                        ::   zzz, zminsal, zmaxsal                 ! local real 
      REAL(wp), DIMENSION(jpvar+1)    ::   zmax
      REAL(wp), DIMENSION(jptst)      ::   zmaxlocal
      LOGICAL                         ::   ll_wrtstp, ll_colruns, ll_wrtruns
      LOGICAL, DIMENSION(jpi,jpj,jpk) ::   llmsk
      CHARACTER(len=20)               ::   clname
      !!----------------------------------------------------------------------
      IF( nstop > 0 .AND. ngrdstop > -1 )   RETURN   !   stpctl was already called by a child grid
      !
      IF( ln_timing )   CALL timing_start( 'stp_ctl' )
      !
      ll_wrtstp  = ( MOD( kt-nit000, sn_cfctl%ptimincr ) == 0 ) .OR. ( kt == nitend )
      ll_colruns = sn_cfctl%l_runstat .AND. ll_wrtstp .AND. jpnij > 1
      ll_wrtruns = sn_cfctl%l_runstat .AND. ll_wrtstp .AND. lwm
      !
      IF( kt == nit000 ) THEN
         !
         IF( lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'stp_ctl : time-stepping control'
            WRITE(numout,*) '~~~~~~~'
         ENDIF
         !                                ! open time.step    ascii file, done only by 1st subdomain
#if defined key_USE_PDAF
         IF( lwm )   CALL ctl_opn( numstp, 'time.step_'//trim(task_str), 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
#else
         IF( lwm )   CALL ctl_opn( numstp, 'time.step', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
#endif
         !
         IF( ll_wrtruns ) THEN
            !                             ! open run.stat     ascii file, done only by 1st subdomain
            CALL ctl_opn( numrun, 'run.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
            !                             ! open run.stat.nc netcdf file, done only by 1st subdomain
            clname = 'run.stat.nc'
            IF( .NOT. Agrif_Root() )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
            istatus = NF90_CREATE( TRIM(clname), NF90_CLOBBER, nrunid )
            istatus = NF90_DEF_DIM( nrunid, 'time', NF90_UNLIMITED, idtime )
            istatus = NF90_DEF_VAR( nrunid, 'abs_ssh_max', NF90_DOUBLE, (/ idtime /), nvarid(1) )
            istatus = NF90_DEF_VAR( nrunid,   'abs_u_max', NF90_DOUBLE, (/ idtime /), nvarid(2) )
            istatus = NF90_DEF_VAR( nrunid,   'abs_v_max', NF90_DOUBLE, (/ idtime /), nvarid(3) )
            istatus = NF90_DEF_VAR( nrunid,       's_min', NF90_DOUBLE, (/ idtime /), nvarid(4) )
            istatus = NF90_DEF_VAR( nrunid,       's_max', NF90_DOUBLE, (/ idtime /), nvarid(5) )
            istatus = NF90_DEF_VAR( nrunid,       't_min', NF90_DOUBLE, (/ idtime /), nvarid(6) )
            istatus = NF90_DEF_VAR( nrunid,       't_max', NF90_DOUBLE, (/ idtime /), nvarid(7) )
            IF( ln_zad_Aimp ) THEN
               istatus = NF90_DEF_VAR( nrunid,   'Cf_max', NF90_DOUBLE, (/ idtime /), nvarid(8) )
               istatus = NF90_DEF_VAR( nrunid,'abs_wi_max',NF90_DOUBLE, (/ idtime /), nvarid(9) )
            ENDIF
            istatus = NF90_ENDDEF(nrunid)
         ENDIF
         !
         ! Set no wet-point indicators (only true if there are no wet-points in
         ! inner domain on respective grids. Note these are surface checks only
         ! - domains containing only cavities are unlikely.) 
         !
         llmsk(     1:nn_hls,:,:) = .FALSE.                                          ! exclude halos from the checked region
         llmsk(Nie0+1:   jpi,:,:) = .FALSE.
         llmsk(:,     1:nn_hls,:) = .FALSE.
         llmsk(:,Nje0+1:   jpj,:) = .FALSE.
         !
         llmsk(Nis0:Nie0,Njs0:Nje0,1) = ssmask(Nis0:Nie0,Njs0:Nje0) == 1._wp         ! define only the inner domain
         l0oce_T = .NOT. ANY( llmsk(:,:,1) )                                         ! no ocean point in the inner domain?
         !
         llmsk(Nis0:Nie0,Njs0:Nje0,1) = umask(Nis0:Nie0,Njs0:Nje0,1) == 1._wp        ! define only the inner domain
         l0oce_U = .NOT. ANY( llmsk(:,:,1) )                                         ! no ocean U-point in the inner domain?
         !
         llmsk(Nis0:Nie0,Njs0:Nje0,1) = vmask(Nis0:Nie0,Njs0:Nje0,1) == 1._wp        ! define only the inner domain
         l0oce_V = .NOT. ANY( llmsk(:,:,1) )                                         ! no ocean V-point in the inner domain?
         IF( jpnij > 1 ) THEN
           j0oce = 0
           IF( l0oce_T ) j0oce(1) = 1
           IF( l0oce_U ) j0oce(2) = 1
           IF( l0oce_V ) j0oce(3) = 1
           CALL mpp_min( "stpctl", j0oce )                                           ! min over the global domain (only 1 if no-sea anywhere)
         ENDIF
      ENDIF
      !
      !                                   !==              write current time step              ==!
      !                                   !==  done only by 1st subdomain at writting timestep  ==!
      IF( lwm .AND. ll_wrtstp ) THEN
         WRITE ( numstp, '(1x, i8)' )   kt
         REWIND( numstp )
      ENDIF
      !                                   !==            test of local extrema           ==!
      !                                   !==  done by all processes at every time step  ==!
      !
      llmsk(     1:nn_hls,:,:) = .FALSE.                                          ! exclude halos from the checked region
      llmsk(Nie0+1:   jpi,:,:) = .FALSE.
      llmsk(:,     1:nn_hls,:) = .FALSE.
      llmsk(:,Nje0+1:   jpj,:) = .FALSE.
      !
      llmsk(Nis0:Nie0,Njs0:Nje0,1) = ssmask(Nis0:Nie0,Njs0:Nje0) == 1._wp         ! define only the inner domain
      !
      IF( ll_wd ) THEN
         zmax(1) = MAXVAL( ABS( ssh(:,:,Kmm) + ssh_ref ), mask = llmsk(:,:,1) )   ! ssh max
      ELSE
         zmax(1) = MAXVAL( ABS( ssh(:,:,Kmm)           ), mask = llmsk(:,:,1) )   ! ssh max
      ENDIF
      llmsk(Nis0:Nie0,Njs0:Nje0,:) = umask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
      zmax(2) = MAXVAL(  ABS( uu(:,:,:,Kmm) ), mask = llmsk )                     ! velocity max (zonal only)
      llmsk(Nis0:Nie0,Njs0:Nje0,:) = vmask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
      zmax(3) = MAXVAL(  ABS( vv(:,:,:,Kmm) ), mask = llmsk )                     ! velocity max (meridional only)
      llmsk(Nis0:Nie0,Njs0:Nje0,:) = tmask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
      zmax(4) = MAXVAL( -ts(:,:,:,jp_sal,Kmm), mask = llmsk )                     ! minus salinity max
      zmax(5) = MAXVAL(  ts(:,:,:,jp_sal,Kmm), mask = llmsk )                     !       salinity max
      IF( ll_colruns .OR. jpnij == 1 ) THEN     ! following variables are used only in the netcdf file
         zmax(6) = MAXVAL( -ts(:,:,:,jp_tem,Kmm), mask = llmsk )                  ! minus temperature max
         zmax(7) = MAXVAL(  ts(:,:,:,jp_tem,Kmm), mask = llmsk )                  !       temperature max
         IF( ln_zad_Aimp ) THEN
            zmax(8) = MAXVAL(   Cu_adv(:,:,:)   , mask = llmsk )                  ! partitioning coeff. max
            llmsk(Nis0:Nie0,Njs0:Nje0,:) = wmask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp
            zmax(9) = MAXVAL(  ABS( wi(:,:,:) ) , mask = llmsk )                  ! implicit vertical vel. max
         ELSE
            zmax(8:9) = 0._wp
         ENDIF
      ELSE
         zmax(6:9) = 0._wp
      ENDIF
      zmax(jpvar+1) = REAL( nstop, wp )                                           ! stop indicator
      !
      !                                   !==               get global extrema             ==!
      !                                   !==  done by all processes if writting run.stat  ==!
      IF( ll_colruns ) THEN
         zmaxlocal(:) = zmax(1:jptst)
         CALL mpp_max( "stpctl", zmax )          ! max over the global domain: ok even if l0oce_[T,U,V] = .true. 
         nstop = NINT( zmax(jpvar+1) )           ! update nstop indicator (now shared among all local domains)
         ! if no ocean point: MAXVAL returns -HUGE => we must overwrite this value to avoid error handling below.
         ! no T-points, globally is unlikely but no U- or V- points can occur with uni-directional channels
         IF( j0oce(2) == 1 )   zmax(2) = 0._wp                                      ! default "valid" values... (in case single-width N-S channel)
         IF( j0oce(3) == 1 )   zmax(3) = 0._wp                                      ! default "valid" values... (in case single-width E-W channel)
      ELSE
         ! Local checks only
         ! if no ocean point: MAXVAL returns -HUGE => we must overwrite this value to avoid error handling below.
         IF( l0oce_T )   zmax(1:jptst) = (/ 0._wp, 0._wp, 0._wp, -1._wp, 1._wp /)   ! default "valid" values...
         IF( l0oce_U )   zmax(2) = 0._wp                                            ! default "valid" values... (in case single-width N-S channel)
         IF( l0oce_V )   zmax(3) = 0._wp                                            ! default "valid" values... (in case single-width E-W channel)
      ENDIF
      !
      zmax(4) = -zmax(4)                              ! move back from max(-zz) to min(zz) : easier to manage! 
      zmax(6) = -zmax(6)                              ! move back from max(-zz) to min(zz) : easier to manage!
      IF( ll_colruns ) zmaxlocal(4) = -zmaxlocal(4)   ! move back from max(-zz) to min(zz) : easier to manage!
      !
      !                                   !==              write "run.stat" files              ==!
      !                                   !==  done only by 1st subdomain at writting timestep  ==!
      IF( ll_wrtruns ) THEN
         WRITE(numrun,9500) kt, zmax(1:jptst)
         IF( jpnij == 1 ) CALL FLUSH(numrun)
         DO ji = 1, jpvar - 2 * COUNT( .NOT. (/ln_zad_Aimp/) )
            istatus = NF90_PUT_VAR( nrunid, nvarid(ji), (/zmax(ji)/), (/kt/), (/1/) )
         END DO
         IF( kt == nitend )   istatus = NF90_CLOSE(nrunid)
      END IF
      !                                   !==               error handling               ==!
      !                                   !==  done by all processes at every time step  ==!
      !
      IF ( ln_SEOS.AND.(rn_b0==0._wp) ) THEN             ! Discard checks on salinity
         zmaxsal =  HUGE(1._wp)                               ! if not used in eos
         zminsal = -HUGE(1._wp)
      ELSE
         zmaxsal = 100._wp
         zminsal =   0._wp
      ENDIF 
      ! 
      IF(  zmax(1) >   20._wp .OR.                  &         ! too large sea surface height ( > 20 m )
         & zmax(2) >   10._wp .OR.                  &         ! too large velocity ( > 10 m/s)
         & zmax(3) >   10._wp .OR.                  &         ! too large velocity ( > 10 m/s)
         & zmax(4) <= zminsal .OR.                  &         ! negative or zero sea surface salinity
         & zmax(5) >= zmaxsal .OR.                  &         ! too large sea surface salinity ( > 100 )
         & zmax(5) <  zminsal .OR.                  &         ! too large sea surface salinity (keep this line for sea-ice)
         & ieee_is_nan( SUM(zmax(1:jptst)) ) .OR.   &         ! NaN encounter in the tests
         & ABS(   SUM(zmax(1:jptst)) ) > HUGE(1._wp) ) THEN   ! Infinity encounter in the tests
         !
         iloc(:,:) = 0
         IF( ll_colruns ) THEN   ! zmax is global, so it is the same on all subdomains -> no dead lock with mpp_maxloc
            ! first: close the netcdf file, so we can read it
            IF( lwm .AND. kt /= nitend )   istatus = NF90_CLOSE(nrunid)
            ! get global loc on the min/max
            llmsk(Nis0:Nie0,Njs0:Nje0,1) = ssmask(Nis0:Nie0,Njs0:Nje0 ) == 1._wp         ! define only the inner domain
            CALL mpp_maxloc( 'stpctl', ABS(ssh(:,:,         Kmm)), llmsk(:,:,1), zzz, iloc(1:2,1) )   ! mpp_maxloc ok if mask = F 
            llmsk(Nis0:Nie0,Njs0:Nje0,:) = umask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
            CALL mpp_maxloc( 'stpctl', ABS( uu(:,:,:,       Kmm)), llmsk(:,:,:), zzz, iloc(1:3,2) )
            llmsk(Nis0:Nie0,Njs0:Nje0,:) = vmask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
            CALL mpp_maxloc( 'stpctl', ABS( vv(:,:,:,       Kmm)), llmsk(:,:,:), zzz, iloc(1:3,3) )
            llmsk(Nis0:Nie0,Njs0:Nje0,:) = tmask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
            CALL mpp_minloc( 'stpctl',      ts(:,:,:,jp_sal,Kmm) , llmsk(:,:,:), zzz, iloc(1:3,4) )
            CALL mpp_maxloc( 'stpctl',      ts(:,:,:,jp_sal,Kmm) , llmsk(:,:,:), zzz, iloc(1:3,5) )
            ! find which subdomain has the max.
            iareamin(:) = jpnij+1   ;   iareamax(:) = 0   ;   iareasum(:) = 0
            DO ji = 1, jptst
               IF( zmaxlocal(ji) == zmax(ji) ) THEN
                  iareamin(ji) = narea   ;   iareamax(ji) = narea   ;   iareasum(ji) = 1
               ENDIF
            END DO
            CALL mpp_min( "stpctl", iareamin )         ! min over the global domain
            CALL mpp_max( "stpctl", iareamax )         ! max over the global domain
            CALL mpp_sum( "stpctl", iareasum )         ! sum over the global domain
         ELSE                    ! find local min and max locations:
            ! if we are here, this means that the subdomain contains some oce points -> no need to test the mask used in maxloc
            llmsk(Nis0:Nie0,Njs0:Nje0,1) = ssmask(Nis0:Nie0,Njs0:Nje0 ) == 1._wp        ! define only the inner domain
            iloc(1:2,1) = MAXLOC( ABS( ssh(:,:,         Kmm)), mask = llmsk(:,:,1) )
            llmsk(Nis0:Nie0,Njs0:Nje0,:) = umask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
            iloc(1:3,2) = MAXLOC( ABS(  uu(:,:,:,       Kmm)), mask = llmsk(:,:,:) )
            llmsk(Nis0:Nie0,Njs0:Nje0,:) = vmask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
            iloc(1:3,3) = MAXLOC( ABS(  vv(:,:,:,       Kmm)), mask = llmsk(:,:,:) )
            llmsk(Nis0:Nie0,Njs0:Nje0,:) = tmask(Nis0:Nie0,Njs0:Nje0,:) == 1._wp        ! define only the inner domain
            iloc(1:3,4) = MINLOC(       ts(:,:,:,jp_sal,Kmm) , mask = llmsk(:,:,:) )
            iloc(1:3,5) = MAXLOC(       ts(:,:,:,jp_sal,Kmm) , mask = llmsk(:,:,:) )
            DO ji = 1, jptst   ! local domain indices ==> global domain indices, excluding halos
               iloc(1:2,ji) = (/ mig(iloc(1,ji),0), mjg(iloc(2,ji),0) /)
            END DO
            iareamin(:) = narea   ;   iareamax(:) = narea   ;   iareasum(:) = 1         ! this is local information
         ENDIF
         !
         WRITE(ctmp1,*) ' stp_ctl: |ssh| > 20 m  or  |U| > 10 m/s  or  S <= 0  or  S >= 100  or  NaN encounter in the tests'
         CALL wrt_line( ctmp2, kt, '|ssh| max', zmax(1), iloc(:,1), iareasum(1), iareamin(1), iareamax(1) )
         CALL wrt_line( ctmp3, kt, '|U|   max', zmax(2), iloc(:,2), iareasum(2), iareamin(2), iareamax(2) )
         CALL wrt_line( ctmp4, kt, '|V|   max', zmax(3), iloc(:,3), iareasum(3), iareamin(3), iareamax(3) )
         CALL wrt_line( ctmp5, kt, 'Sal   min', zmax(4), iloc(:,4), iareasum(4), iareamin(4), iareamax(4) )
         CALL wrt_line( ctmp6, kt, 'Sal   max', zmax(5), iloc(:,5), iareasum(5), iareamin(5), iareamax(5) )
         IF( Agrif_Root() ) THEN
            WRITE(ctmp7,*) '      ===> output of last computed fields in output.abort* files'
         ELSE
            WRITE(ctmp7,*) '      ===> output of last computed fields in '//TRIM(Agrif_CFixed())//'_output.abort* files'
         ENDIF
         !
         CALL dia_wri_state( Kmm, 'output.abort' )     ! create an output.abort file
         !
         IF( ll_colruns .OR. jpnij == 1 ) THEN   ! all processes synchronized -> use lwp to print in opened ocean.output files
            IF(lwp .AND. ln_zad_Aimp) WRITE(numout,*) 'max partitioning coefficient and max wi: ',zmax(8), zmax(9)
            IF( .NOT. ll_colruns .AND. lwm )   istatus = NF90_CLOSE(nrunid)
            IF(lwp) THEN   ;   CALL ctl_stop( ctmp1, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ' ', ctmp7 )
            ELSE           ;   nstop = MAX(1, nstop)   ! make sure nstop > 0 (automatically done when calling ctl_stop)
            ENDIF
         ELSE                                    ! only mpi subdomains with errors are here -> STOP now
            IF( lwm )   istatus = NF90_CLOSE(nrunid)
            IF(ln_zad_Aimp) WRITE(*,*) 'max coeff and max |wi|: ',narea, zmax(8), zmax(9)
            CALL ctl_stop( 'STOP', ctmp1, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ' ', ctmp7 )
         ENDIF
         !
      ENDIF
      !
      IF( nstop > 0 ) THEN                                                  ! an error was detected and we did not abort yet...
         ngrdstop = Agrif_Fixed()                                           ! store which grid got this error
         IF( .NOT. ll_colruns .AND. jpnij > 1 )   CALL ctl_stop( 'STOP' )   ! we must abort here to avoid MPI deadlock
      ENDIF
      !
9500  FORMAT(' it :', i8, '    |ssh|_max: ', D23.16, ' |U|_max: ', D23.16, ' |V|_max: ', D23.16, ' S_min: ', D23.16,' S_max: ', D23.16)
      !
      IF( ln_timing )   CALL timing_stop( 'stp_ctl' )
      !
   END SUBROUTINE stp_ctl


   SUBROUTINE wrt_line( cdline, kt, cdprefix, pval, kloc, ksum, kmin, kmax )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE wrt_line  ***
      !!
      !! ** Purpose :   write information line
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*),      INTENT(  out) ::   cdline
      CHARACTER(len=*),      INTENT(in   ) ::   cdprefix
      REAL(wp),              INTENT(in   ) ::   pval
      INTEGER, DIMENSION(3), INTENT(in   ) ::   kloc
      INTEGER,               INTENT(in   ) ::   kt, ksum, kmin, kmax
      !
      CHARACTER(len=80) ::   clsuff
      CHARACTER(len=9 ) ::   clkt, clsum, clmin, clmax
      CHARACTER(len=9 ) ::   cli, clj, clk
      CHARACTER(len=1 ) ::   clfmt
      CHARACTER(len=4 ) ::   cl4   ! needed to be able to compile with Agrif, I don't know why
      INTEGER           ::   ifmtk
      !!----------------------------------------------------------------------
      WRITE(clkt , '(i9)') kt
      
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpnij  ,wp))) + 1     ! how many digits to we need to write ? (we decide max = 9)
      !!! WRITE(clsum, '(i'//clfmt//')') ksum                   ! this is creating a compilation error with AGRIF
      cl4 = '(i'//clfmt//')'   ;   WRITE(clsum, cl4) ksum
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(MAX(1,jpnij-1),wp))) + 1    ! how many digits to we need to write ? (we decide max = 9)
      cl4 = '(i'//clfmt//')'   ;   WRITE(clmin, cl4) kmin-1
                                   WRITE(clmax, cl4) kmax-1
      !
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpiglo,wp))) + 1      ! how many digits to we need to write jpiglo? (we decide max = 9)
      cl4 = '(i'//clfmt//')'   ;   WRITE(cli, cl4) kloc(1)      ! this is ok with AGRIF
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpjglo,wp))) + 1      ! how many digits to we need to write jpjglo? (we decide max = 9)
      cl4 = '(i'//clfmt//')'   ;   WRITE(clj, cl4) kloc(2)      ! this is ok with AGRIF
      !
      IF( ksum == 1 ) THEN   ;   WRITE(clsuff,9100) TRIM(clmin)
      ELSE                   ;   WRITE(clsuff,9200) TRIM(clsum), TRIM(clmin), TRIM(clmax)
      ENDIF
      IF(kloc(3) == 0) THEN
         ifmtk = INT(LOG10(REAL(jpk,wp))) + 1                   ! how many digits to we need to write jpk? (we decide max = 9)
         clk = REPEAT(' ', ifmtk)                               ! create the equivalent in blank string
         WRITE(cdline,9300) TRIM(ADJUSTL(clkt)), TRIM(ADJUSTL(cdprefix)), pval, TRIM(cli), TRIM(clj), clk(1:ifmtk), TRIM(clsuff)
      ELSE
         WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpk,wp))) + 1      ! how many digits to we need to write jpk? (we decide max = 9)
         !!! WRITE(clk, '(i'//clfmt//')') kloc(3)               ! this is creating a compilation error with AGRIF
         cl4 = '(i'//clfmt//')'   ;   WRITE(clk, cl4) kloc(3)   ! this is ok with AGRIF
         WRITE(cdline,9400) TRIM(ADJUSTL(clkt)), TRIM(ADJUSTL(cdprefix)), pval, TRIM(cli), TRIM(clj),    TRIM(clk), TRIM(clsuff)
      ENDIF
      !
9100  FORMAT('MPI rank ', a)
9200  FORMAT('found in ', a, ' MPI tasks, spread out among ranks ', a, ' to ', a)
9300  FORMAT('kt ', a, ' ', a, ' ', 1pg11.4, ' at i j   ', a, ' ', a, ' ', a, ' ', a)
9400  FORMAT('kt ', a, ' ', a, ' ', 1pg11.4, ' at i j k ', a, ' ', a, ' ', a, ' ', a)
      !
   END SUBROUTINE wrt_line


   !!======================================================================
END MODULE stpctl
