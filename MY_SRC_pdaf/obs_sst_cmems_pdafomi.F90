!$Id: obs_TYPE_pdafomi_TEMPLATE.F90 579 2020-11-23 07:32:00Z lnerger $
!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!!
!! Observation type: Regional SST observations from CMEMS
!!
!! The subroutines in this module are for the particular handling of
!! satellite SST observations in regional configurations that use a
!! regular lon/lat model grid. This type of the model grid is used for
!! the efficient computation of the indices of grid points surrounding
!! an observation location. The observation module can handle the case
!! that model and observation grid only overlap partially. Both
!! interpolation to the observation grid and superobbing by averaging
!! observations on the model grid are supported.
!!
!! A particularity of the observation type used here is that it is stored
!! in short integer form. For this integer values are read from the file
!! and converted into real.
!!
!! The routines are called by the different call-back routines of PDAF.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
module obs_sst_cmems_pdafomi

  USE mod_kind_pdaf
  use parallel_pdaf, &
       only: mype_filter    ! Rank of filter process
  use PDAFomi, &
       only: obs_f, obs_l   ! Declaration of observation data types
 
  implicit none
  save

  ! Variables which are inputs to the module (usually set in init_pdaf)
  logical :: assim_sst_cmems = .false.  !< Whether to assimilate this data type
  real(pwp) :: rms_obs_sst_cmems = 0.8  !< Observation error standard deviation (for constant errors)
  real(pwp) :: lradius_sst_cmems = 1.0  !< Localization cut-off radius
  real(pwp) :: sradius_sst_cmems = 1.0  !< Support radius for weight function
  integer :: mode_sst_cmems = 0         !< Observation mode: 
                                        !< (0) linear interpolation
                                        !< (1) super-obbing: average 4 observation values
  integer :: time_sst_cmems = 0         !< Time at which the daily observations are assimilated
                                        !< (0) mightnight, (12) noon
  character(len=3) :: dist_sst_cmems = 'geo'  ! Type of distance computation: 
                                        !< (gp) for Cartesian distance in unit of grid points
                                        !< (geo) for geographic distance in km
                                        !<  Note: The implementation assumes a regular lat/lon grid
  character(lc) :: path_sst_cmems = '.' !< Path to observation file
  character(lc) :: file_sst_cmems = ''  !< Filename for observations
  character(lc) :: varname_sst_cmems = 'adjusted_sea_surface_temperature'  !< Name of observation variable in file

! *** Variables used inside the module
  integer, private :: observation_mode        !< Observation mode (alias for mode_sst_cmems)
  character(len=3) :: dist_obs                !< Type of distance computation (alias for dist_sst_cmems)
  character(lc) :: obsname = 'OBS_SST_CMEMS'  !< Observation name string

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  type(obs_f), target, public :: thisobs      ! full observation
  type(obs_l), target, public :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

contains

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  subroutine init_dim_obs_sst_cmems(step, dim_obs)

    use netcdf
    use PDAFomi, &
         only: PDAFomi_gather_obs, PDAFomi_get_interp_coeff_lin
    use assimilation_pdaf, &
         only: filtertype, screen
    use statevector_pdaf, &
         only: id, sfields
    use parallel_pdaf, &
         only: mype_filter, npes_filter
    use io_pdaf, &
         only: check
    use nemo_pdaf, &
         only: lat1_p, lon1_p, nlats=>nj_p, nlons=>ni_p, &
         idx_nwet, use_wet_state, nlei, nlej, calc_date, deg2rad

    implicit none

! *** Arguments ***
    integer, intent(in)    :: step       !< Current time step
    integer, intent(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    logical :: debug = .false.               ! Activate debugging output for index calculations
    logical :: obsgrid_p = .true.            ! Flag whether the observation grid includes PE-local sub-domain
    integer :: i, j, cnt                     ! Counters
    integer :: ido_start, ido_end            ! Counters
    integer :: idm_start, idm_end            ! Counters
    integer :: dim_obs_p                     ! Number of process-local observations
    real(pwp), allocatable :: obs_p(:)       ! PE-local observation vector
    real(pwp), allocatable :: ivar_obs_p(:)  ! PE-local inverse observation error variance
    real(pwp), allocatable :: ocoord_p(:,:)  ! PE-local observation coordinates 
    logical :: doassim_now=.true.            ! Whether we assimilate the observation at the current time
    integer(4) :: status                     ! Status flag for availability of observations
    character(len=100) :: file_full          ! filename including path
    character(len=2) :: strday               ! day as string
    integer(4) :: ncid, dimid, lonid, latid, varid       ! nc file IDs
    integer(4) :: startv(3), cntv(3)                     ! Index arrays for reading from nc file
    integer(4) :: dim_olat, dim_olon                     ! Grid dimensions read from file
    integer(4), allocatable :: obs_from_file(:,:)        ! observation field read from file 
    real(pwp), allocatable :: lon_obs(:), lat_obs(:)     ! Obs. coordinates read from file
    real(pwp), allocatable :: lon_model(:), lat_model(:) ! Longitude/latitude of model in radians
    integer :: iderr                         ! Error flag for determining indices
    integer(4) :: ido_n, ido_e, ido_s, ido_w ! Obs. ID limits N/E/S/W for model grid
    integer :: idm_n, idm_e, idm_s, idm_w    ! Model ID limits NESW 
    real(pwp) :: wlonM, elonM, nlatM, slatM  ! Coordinate limits of model grid
    real(pwp) :: dlonM, dlatM                ! Model grid spacing
    real(pwp) :: dlonO, dlatO                ! Observation grid spacing
    real(pwp) :: latM_limit                  ! Comparison limit in latitude for model coordinate
    real(pwp) :: lonM, latM                  ! Longitude/latitude of a model grid point
    real(pwp) :: gcoords(4,2)                ! Grid point coordinates for computing interpolation coeffs
    integer(4) :: obsflag                    ! Count observation in direct vicinity
    integer(4) :: cntobs(4)                  ! Count grid points with 0 to 4 obs. neighbours
    integer(4) :: obs_sum                    ! Sum of observation integer values
    integer :: sgn_olat                      ! Orientation of latitude Obs.: -1 for north-south/+1 for south-north
    integer :: sgn_mlat                      ! Orientation of latitude Model: -1 for north-south/+1 for south-north
    real(pwp) :: rdate                       ! Current date
    integer(4) :: year, month, iday          ! Current year, month, day (iday is step read from observation file) 
    integer :: id_obs                        ! Index of observation field in state vector
    character(lc) :: varname_lon             ! Name of longitude coordinate variable in file
    character(lc) :: varname_lat             ! Name of latitude coordinate variable in file
    character(lc) :: varname_obs             ! Name of observation variable in file
    real(pwp) :: rms_obs                     ! Obs. error standard deviation
    integer :: missing_value                 ! Missing value above which observations are valid
    character(len=2) :: region               ! Region for which the data is used ('no', 'ba', 'nb')
    real(pwp) :: limcoords(3)                ! Limiting coordinates according to region
    real(pwp), parameter :: sst_scale = 0.01 ! Scaling factor to convert file value to degC


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    if (mype_filter==0) &
         write (*,'(a,4x,a,a)') 'NEMO-PDAF', 'Assimilate observations - ', trim(obsname)

    ! Specify region and limiting coordinates
    region = 'nb'                     ! 'nb'= North and Baltic Seas - no exclusion
    limcoords(1) = 57.0 * deg2rad    ! north/south limit in Skagerrak
    limcoords(2) = 9.4 * deg2rad      ! east/west limit over Denmark; use outh of limcoords(1) for 'no'
    limcoords(3) = 15.0 * deg2rad     ! east/west limit over Sweden; use north of limcoords(1) for 'ba'

    ! Initialize generic variables (used to keep codes generic)
    id_obs = id%temp                         ! Index of observation field in state vector
    dist_obs = dist_sst_cmems                ! Type of distance computation
    observation_mode = mode_sst_cmems        ! Whether to use interpolation or super-obbing
    varname_lon = 'lon'                      ! Name of longitude variable in file
    varname_lat = 'lat'                      ! Name of latitude variable in file
    varname_obs = varname_sst_cmems          ! Name of observation variable in file
    rms_obs = rms_obs_sst_cmems              ! Obs. error standard deviation
    missing_value = -10000                   ! Missing value in observation file
    file_full = trim(path_sst_cmems)//trim(file_sst_cmems)   ! File name including path

    ! Store whether to assimilate this observation type (used in routines below)
    if (assim_sst_cmems) thisobs%doassim = 1

    ! Specify type of distance computation
    if (trim(dist_obs) == 'gp') then
       if (mype_filter==0) write (*,'(a,4x,a)') 'NEMO-PDAF', '--- use Cartesian grid point distances'
       thisobs%disttype = 0   ! 0=Cartesian
    else
       if (mype_filter==0) write (*,'(a,4x,a)') 'NEMO-PDAF', '--- use geographic distances'
       thisobs%disttype = 2   ! 2=Geographic
    end if

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! In case of MPI parallelization restrict observations to sub-domains
    if (npes_filter>1) thisobs%use_global_obs = 0


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Determine current day
    call calc_date(step-1, rdate)
    year = floor(rdate/10000.0_pwp)
    month = floor((rdate-real(year*10000))/100.0_pwp)
    iday = floor(rdate-real(year*10000)-real(month*100))

    doassim: if (doassim_now) then

       ! Only execute this if we assimilate observations at this hour

       ! The SST data can can be downloaded as follows:
       ! python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu 
       !    --service-id SST_EUR_SST_L3S_NRT_OBSERVATIONS_010_009_a-TDS 
       !    --product-id METEOFRANCE-EUR-SST_L3MULTISENSOR_NRT-OBS_FULL_TIME_SERIE 
       !    --longitude-min -4.5 --longitude-max 30.5 --latitude-min 48.5 --latitude-max 66 
       !    --date-min "2018-10-01 00:00:00" --date-max "2018-10-31 00:00:00" 
       !    --variable adjusted_sea_surface_temperature 
       !    --out-name sst_multi_201810.nc --user <USERNAME> --pwd <PASSWD>

       ! read observation values and their coordinates

       if (mype_filter==0) then 
          write (*,'(a, 4x,a,i3,a)') 'NEMO-PDAF', 'Read observations for day ', &
               iday,' from file:'
          write (*,'(a, 4x,a)') 'NEMO-PDAF', trim(file_full)
          write (*, '(a,4x,a,a)') 'NEMO-PDAF', '--- name of observation file variable: ', trim(varname_obs)
       end if

       ! Open the file. NF90_NOWRITE tells netCDF to have read-only access to file.
       call check( nf90_open(file_full, NF90_NOWRITE, ncid) )

       ! Read dimensions of observation grid
       call check( nf90_inq_dimid(ncid, trim(varname_lat), dimid) )
       call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olat) )
       call check( nf90_inq_dimid(ncid, trim(varname_lon), dimid) )
       call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olon) )

       ! Allocate arrays
       allocate(obs_from_file(dim_olon, dim_olat))
       allocate(lon_obs(dim_olon), lat_obs(dim_olat))

       ! Get variable IDs and read data
       call check( nf90_inq_varid(ncid, trim(varname_obs), varid) )
       call check( nf90_inq_varid(ncid, trim(varname_lon), lonid) )
       call check( nf90_inq_varid(ncid, trim(varname_lat), latid) )

       call check( nf90_get_var(ncid, lonid, lon_obs) )
       call check( nf90_get_var(ncid, latid, lat_obs) )

       ! Read observation values
       ! They are in deg C but have to be scaled by sst_scale (1/100).
       startv(1) = 1 ! lon
       startv(2) = 1 ! lat
       startv(3) = iday ! time
       cntv(1) = dim_olon
       cntv(2) = dim_olat
       cntv(3) = 1
       call check( nf90_get_var(ncid, varid, obs_from_file, start=startv, count=cntv) )

       ! Close the file
       call check( nf90_close(ncid) )

       ! Note: The implementated coordinate handling to find grid points
       ! neighboring an observation assumes a regular lat/lon grid
  
       ! *** Convert observation coordinates to radians ***
       lon_obs = lon_obs * deg2rad
       lat_obs = lat_obs * deg2rad

       ! *** Store model coordinates in radians ***
       allocate(lon_model(nlons))
       allocate(lat_model(nlats))
       lon_model = lon1_p * deg2rad
       lat_model = lat1_p * deg2rad

       cartdist: if (trim(dist_obs) == 'gp') then

          ! Set grid point coordinates for Cartesian distance computation

          dlatM = lat_model(2) - lat_model(1)  ! Model grid spacing in latitude
          dlonM = lon_model(2) - lon_model(1)  ! Model grid spacing in longitude

          do i = 1, dim_olon
             lon_obs(i) = (lon_obs(i) - lon_model(1))/dlonM + 1.0
          end do
          do j = 1, dim_olat
             lat_obs(j) = (lat_obs(j) - lat_model(1))/dlatM + 1.0
          end do

          do i = 1, nlons
             lon_model(i) = i + nlei - 1
          end do
          do j = 1, nlats
             lat_model(j) = j + nlej - 1
          end do

       end if cartdist


! ****************************************
! *** Exclude data for specific region ***
! ****************************************

       if (region=='ba') then
          ! Baltic Sea (exclude Skagerrak/Kattegat and all west of Denmark)

          if (mype_filter==0) &
               write (*,'(8x,a)') '--- Exclude observations in North Sea'

          do j = 1, dim_olat
             do i = 1, dim_olon
                if (lon_obs(i)<limcoords(2) .or. (lon_obs(i)<limcoords(3) &
                     .and. lat_obs(j)>limcoords(1))) then
                   obs_from_file(i,j) = missing_value
                end if
             end do
          end do
       elseif (region=='no') then
          ! North Sea (exclude Baltic except Skagerrak/Kattegat north of limcoords(1))

          if (mype_filter==0) &
               write (*,'(8x,a)') '--- Exclude observations in Baltic Sea'

          do j = 1, dim_olat
             do i = 1, dim_olon
                if ((lon_obs(i)>=limcoords(2) .and. lat_obs(j)<=limcoords(1)) &
                     .or. lon_obs(i)>limcoords(3)) then
                   obs_from_file(i,j) = missing_value
                end if
             end do
          end do
       end if


! *******************************************************************
! *** Determine indices of observation grid inside the model grid ***
! *******************************************************************

       ! Set index limits for model grid
       idm_n = nlats
       idm_s = 1
       idm_w = 1
       idm_e = nlons

       ! boundaries and spacing for model grid
       nlatM = lat_model(idm_n)             ! Northern latitude limit of the grid box
       slatM = lat_model(idm_s)             ! Southern latitude limit of the grid box
       wlonM = lon_model(idm_w)             ! Western longitude limit of the grid box 
       elonM = lon_model(idm_e)             ! Eastern longitude limit of the grid box
       dlatM = lat_model(2) - lat_model(1)  ! Model grid spacing in latitude
       dlonM = lon_model(2) - lon_model(1)  ! Model grid spacing in longitude
       sgn_mlat = int(sign(1.0,dlatM))      ! Orientation of latitudinal direction (-1: north-to-south)

       ! Observation grid spacing and orientation
       dlatO = lat_obs(2) - lat_obs(1)      ! Observation grid spacing in latitude
       dlonO = lon_obs(2) - lon_obs(1)      ! Observation grid spacing in longitude
       sgn_olat = int(sign(1.0,dlatO))      ! Orientation of latitudinal direction (-1: north-to-south)

       if (debug) then
          write (*,*) 'dim_olat/lon', dim_olat, dim_olon
          write (*, *) 'idm all', idm_w, idm_e, idm_n, idm_s
          write (*,'(a,4f12.5)') 'mcoords limits',wlonM, elonM, nlatM, slatM
          if (sgn_mlat<0) then
             write (*,'(a,4f12.5)') 'ocoords limits', lon_obs(1), lon_obs(dim_olon), lat_obs(1), lat_obs(dim_olat)
          else
             write (*,'(a,4f12.5)') 'ocoords limits', lon_obs(1), lon_obs(dim_olon), lat_obs(dim_olat), lat_obs(1)
          end if
          write (*,*) 'sgn_olat/sgn_mlat', sgn_olat, sgn_mlat
          write (*,'(a,2es10.2,a,2es10.2)') 'dlat/lon: M:', dlatM, dlonM, ', O:', dlatO, dlonO
       end if


       ! Initialize error flag
       iderr = 0

       ! Compute obs. coordinate indices for model grid limits 
       ! If needed adapt index limits for model grid
       ido_w = ceiling(abs(wlonM - lon_obs(1)) / dlonO)+1

       if (ido_w<2) then
          if (debug) write (*,*) 'NEMO-PDAF ', 'reset ido_w'
          ido_w = 2
          idm_w = ceiling((lon_obs(ido_w) - wlonM) / dlonM) + 1
       elseif (ido_w>dim_olon) then
          iderr = 1
       endif
       if (observation_mode==0 .and. iderr==0) then
          idm_w = ceiling((lon_obs(ido_w) - wlonM) / dlonM) + 1
       end if

       ido_e = ceiling((elonM - lon_obs(1)) / dlonO)
       if (ido_e>=dim_olon) then
          if (debug) write (*,*) 'NEMO-PDAF ', 'reset ido_e'
          ido_e = dim_olon-1
          idm_e = floor((lon_obs(dim_olon) - wlonM) / dlonM) + 1
       elseif (ido_e<1) then
          iderr = 2
       endif
       if (observation_mode==0 .and. iderr==0) then
          idm_e = floor((lon_obs(ido_e) - wlonM) / dlonM) + 1
       end if

       if (sgn_olat > 0) then
          ido_n = ceiling((nlatM - lat_obs(1)) / dlatO)
          if (ido_n>dim_olat) then
             if (debug) write (*,*) 'NEMO-PDAF ', 'reset ido_n'
             ido_n = dim_olat-1
             idm_n = ceiling((lat_obs(ido_n) - nlatM) / dlatM) + 1
          elseif (ido_n<1) then
             iderr = 3
          end if
          if (observation_mode==0 .and. iderr==0) then
             if (sgn_mlat > 0) then
                idm_n = floor(abs(lat_obs(ido_n) - slatM) / dlatM) + 1
             else
                idm_n = ceiling((lat_obs(ido_n) - nlatM) / dlatM) + 1
             end if
          end if

          ido_s = ceiling(abs(slatM - lat_obs(1)) / dlatO) + 1
          if (ido_s<2) then
             if (debug) write (*,*) 'NEMO-PDAF ', 'reset ido_s'
             ido_s = 2
             if (sgn_mlat > 0) then
                idm_s = floor(abs(lat_obs(ido_s) - slatM) / dlatM) + 1
             else
                idm_s = floor(abs(lat_obs(ido_s) - nlatM) / dlatM) + 1
             end if
          elseif (ido_s>dim_olat) then
             iderr = 4
          end if
          if (observation_mode==0 .and. iderr==0) then
             if (sgn_mlat > 0) then
                idm_s = ceiling(abs(lat_obs(ido_s) - slatM) / dlatM) + 1
             else
                idm_s = floor((lat_obs(ido_s) - nlatM) / dlatM) + 1
             end if
          end if
       else
          ido_n = ceiling((nlatM - lat_obs(1)) / dlatO)+1
          if (ido_n<2) then 
             if (debug) write (*,*) 'NEMO-PDAF ', 'reset ido_n'
             ido_n = 2
             idm_n = floor((lat_obs(ido_n) - nlatM) / dlatM) + 1
          elseif (ido_n>dim_olat) then
             iderr = 5
          endif
          if (observation_mode==0 .and. iderr==0) then
             if (sgn_mlat > 0) then
                idm_n = floor((lat_obs(ido_n) - slatM) / dlatM) + 1
             else
                idm_n = ceiling((lat_obs(ido_n) - nlatM) / dlatM) + 1
             end if
          end if

          ido_s = floor((slatM - lat_obs(1)) / dlatO)
          if (ido_s>dim_olat) then
             if (debug) write (*,*) 'NEMO-PDAF ', 'reset ido_s'
             ido_s = dim_olat-1
             idm_s = ceiling((lat_obs(dim_olat) - nlatM) / dlatM)+1
          elseif (ido_s<1) then
             iderr = 6
          endif
          if (observation_mode==0 .and. iderr==0) then
             if (sgn_mlat > 0) then
                idm_s = ceiling(abs(lat_obs(ido_s) - slatM) / dlatM)
             else
                idm_s = floor(abs(lat_obs(ido_s) - nlatM) / dlatM)+1
             end if
          end if
       end if

       if (iderr/=0) then
          ! Observation grid and PE-local sub-domain do not overlap
          obsgrid_p = .false.
       else
          ! Overlapping observation and PE-local grids
          obsgrid_p = .true.
       end if

       if (debug) then
          if (.not. obsgrid_p) then
             write (*,*) 'Observations not overlapping with model grid, case', iderr
          else
             write (*,*) 'Observations do overlap with model grid'
          end if
          write (*,'(a,2x,4i10)')   'ido in WENS', ido_w, ido_e, ido_n, ido_s
          write (*,'(a,4f12.5)') 'ocoords WENS in ', &
               lon_obs(ido_w), lon_obs(ido_e), lat_obs(ido_n), lat_obs(ido_s)
          write (*,'(a,4f12.5)') 'ocoords WENS out', &
               lon_obs(ido_w-1), lon_obs(ido_e+1), lat_obs(ido_n+sgn_olat), lat_obs(ido_s-sgn_olat)

          write (*,'(a,2x,4i10)') 'idm new limits', idm_w, idm_e, idm_n, idm_s
          if (observation_mode==0) &
               write (*,'(a,4f12.5)') 'mcoords out     ', lon_model(idm_w-1), &
               lon_model(idm_e+1), lat_model(idm_n+sgn_mlat), lat_model(idm_s-sgn_mlat)
          write (*,'(a,4f12.5)') 'mcoords in      ', lon_model(idm_w), &
               lon_model(idm_e), lat_model(idm_n), lat_model(idm_s)
          if (sgn_mlat>0) then
             write (*,*) 'limits SN', 1, nlats, 'new', idm_s, idm_n
          else
             write (*,*) 'limits NS', 1, nlats, 'new', idm_n, idm_s
          end if
          write (*,*) 'limits WE', 1, nlons, 'new', idm_w, idm_e
       end if


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

       obsmodeA: if (observation_mode==0 .and. obsgrid_p) then

          ! *** Linear interpolation ***

          if (mype_filter==0) &
               write (6,'(a,4x,a)') 'NEMO-PDAF', '--- use observations with linear interpolation'

          ! *** Count valid observations that lie within the model grid ***
          cnt = 0

          ! Set start and end index for latitude
          if (sgn_olat == 1) then
             ido_start = ido_s
             ido_end = ido_n
          else
             ido_start = ido_n
             ido_end = ido_s
          end if

          ! Set limit index for computing model latitude index
          if (sgn_mlat>0) then
             latM_limit = slatM
          else
             latM_limit = nlatM
          end if

          ! Loop through observation grid
          do j = ido_start, ido_end
             do i = ido_w, ido_e
                if (obs_from_file(i,j) > missing_value) then

                   ! find model grid point indices corresponding to observation point
                   idm_w = floor((lon_obs(i) - wlonM) / dlonM) + 1
                   idm_e = idm_w + 1
                   idm_n = ceiling(abs(lat_obs(j) - latM_limit) / dlatM) + 1
                   idm_s = idm_n - sgn_mlat

                   ! Check whether all of the neighboring grid points are wet - then use the observation
                   if (idx_nwet(idm_w, idm_n)>0 .and. idx_nwet(idm_e, idm_n)>0 .and. &
                        idx_nwet(idm_w, idm_s)>0 .and. idx_nwet(idm_e, idm_s) > 0) then
                      cnt = cnt + 1
                   end if

                endif
             end do
          end do

       elseif (obsgrid_p) then obsmodeA

          ! *** Super-Obbing ***

          if (mype_filter==0) &
               write (6,'(a,4x,a)') 'NEMO-PDAF', '--- use observations with super-obbing'

          ! *** Count valid observations that lie within the grid ***
          cnt = 0
          cntobs = 0

          ! Set start and end index for latitude
          if (sgn_olat == 1) then
             idm_start = idm_s
             idm_end = idm_n
          else
             idm_start = idm_n
             idm_end = idm_s
          end if

          ! Loop through model grid
          do j = idm_start, idm_end
             do i = idm_w, idm_e

                ! Model grid point coordinates
                latM = lat_model(j)
                lonM = lon_model(i)
             
                ! Compute observation grid point indices
                ido_w = ceiling(abs(lonM - lon_obs(1)) / dlonO)
                ido_e = ceiling(abs(lonM - lon_obs(1)) / dlonO)+1
                ido_n = ceiling(abs(latM - lat_obs(1)) / dlatO)+1
                ido_s = ceiling(abs(latM - lat_obs(1)) / dlatO)
             
                ! For wet grid points: Check whether there are valid
                ! observations at the observation grid points around the point
                if (idx_nwet(i, j)>0.0) then

                   if (ido_w<1 .or. ido_w>dim_olon .or. ido_e<1 .or. ido_e>dim_olon &
                        .or. ido_n<1 .or. ido_n>dim_olat .or. ido_s<1 .or. ido_s>dim_olat) &
                        write (*,*) 'ido out of range', ido_w, ido_e, ido_n, ido_s

                   obsflag = 0
                   if (obs_from_file(ido_w,ido_s) > missing_value) obsflag = obsflag + 1
                   if (obs_from_file(ido_e,ido_s) > missing_value) obsflag = obsflag + 1
                   if (obs_from_file(ido_w,ido_n) > missing_value) obsflag = obsflag + 1
                   if (obs_from_file(ido_e,ido_n) > missing_value) obsflag = obsflag + 1

                   ! Count how many of the neighboring observations are valid
                   if (obsflag==4) cntobs(4) = cntobs(4) + 1
                   if (obsflag==3) cntobs(3) = cntobs(3) + 1
                   if (obsflag==2) cntobs(2) = cntobs(2) + 1
                   if (obsflag==1) cntobs(1) = cntobs(1) + 1

                   ! Count grid points with observations in direct vicinity
                   if (obsflag>0) cnt = cnt+1
                end if

             end do
          end do
          write (6, '(a,8x, a,4i8)') 'NEMO-PDAF', 'grid points with 1/2/3/4 neighbour obs.', cntobs(1:4)

       else obsmodeA
          ! Observation and model grids do not overlap
          cnt = 0
       end if obsmodeA

       ! Set observation dimension
       dim_obs_p = cnt
       dim_obs = cnt 

       if (npes_filter==1) then
          write (6,'(a, 4x, a, a, i7)') 'NEMO-PDAF', '--- number of observations from ', trim(obsname), ': ', dim_obs
       else
          if (screen>2) then
             write (6,'(a, 4x, a, i4, 2x, a, a, i7)') 'NEMO-PDAF', 'PE', mype_filter, &
                  '--- number of observations from ', trim(obsname), ': ', dim_obs
          end if
       end if

    else doassim

       ! This is for the case that we do not assimilate this data at the current hour
       dim_obs = 0
       dim_obs_p = 0

    end if doassim


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations                 ***
    ! *** Initialize process local index array                        ***
    if (dim_obs_p > 0) then

       allocate(obs_p(dim_obs_p))
       allocate(ivar_obs_p(dim_obs_p))
       allocate(ocoord_p(thisobs%ncoord, dim_obs_p))

       obsmodeB: if (observation_mode==0) then

          ! *** Linear interpolation ***

          ! Allocate process-local index array
          allocate(thisobs%id_obs_p(4, dim_obs_p))

          ! Allocate array of interpolation coefficients. As ID_OBS_P, the number
          ! of rows corresponds to the number of grid points used in the interpolation
          allocate(thisobs%icoeff_p(4, dim_obs_p))

          cnt = 0

          do j = ido_start, ido_end
             do i = ido_w, ido_e

                ! find grid point indices corresponding to observation point
                idm_w = floor((lon_obs(i) - wlonM) / dlonM) + 1
                idm_e = idm_w + 1
                idm_n = ceiling(abs(lat_obs(j) - latM_limit) / dlatM) + 1
                idm_s = idm_n - sgn_mlat

                haveobs: if (obs_from_file(i,j) > missing_value) then

                   ! find grid point indices corresponding to observation point
                   idm_w = floor((lon_obs(i) - wlonM) / dlonM) + 1
                   idm_e = idm_w + 1
                   idm_n = ceiling(abs(lat_obs(j) - latM_limit) / dlatM) + 1
                   idm_s = idm_n - sgn_mlat

                   ! Use obs. if all of the neighboring grid points are wet
                   wetpoint: if (idx_nwet(idm_w, idm_n)>0 .and. idx_nwet(idm_e, idm_n)>0 .and. &
                        idx_nwet(idm_w, idm_s)>0 .and. idx_nwet(idm_e, idm_s) > 0) then

                      cnt = cnt + 1

                      ! Set indices of 4 grid points neightbors of the observation
                      if (use_wet_state==1 .or. use_wet_state==2) then
                         thisobs%id_obs_p(1, cnt) = idx_nwet(idm_w,idm_s) + sfields(id_obs)%off
                         thisobs%id_obs_p(2, cnt) = idx_nwet(idm_e,idm_s) + sfields(id_obs)%off
                         thisobs%id_obs_p(3, cnt) = idx_nwet(idm_w,idm_n) + sfields(id_obs)%off
                         thisobs%id_obs_p(4, cnt) = idx_nwet(idm_e,idm_n) + sfields(id_obs)%off
                      else
                         thisobs%id_obs_p(1, cnt) = idm_w + nlons*(idm_s-1) + sfields(id_obs)%off
                         thisobs%id_obs_p(2, cnt) = idm_e + nlons*(idm_s-1) + sfields(id_obs)%off
                         thisobs%id_obs_p(3, cnt) = idm_w + nlons*(idm_n-1) + sfields(id_obs)%off
                         thisobs%id_obs_p(4, cnt) = idm_e + nlons*(idm_n-1) + sfields(id_obs)%off
                      end if

                      ! Store observation value and coordinates
                      obs_p(cnt) = real(obs_from_file(i, j)) * sst_scale
                      ocoord_p(1, cnt) = lon_obs(i)
                      ocoord_p(2, cnt) = lat_obs(j)

                      ! *** Determine interpolation coefficients ***

                      ! Determine coordinates of grid points around observation
                      ! Order of coefficients:  (3) ---- (4)          
                      !                          |        |
                      !                         (1) ---- (2)
                      ! Only 4 coordinate values are required for bi-linear interpolation
                      gcoords(1,1) = lon_model(idm_w)
                      gcoords(1,2) = lat_model(idm_s)
                      gcoords(2,1) = lon_model(idm_e)
                      gcoords(3,2) = lat_model(idm_n)

                      ! Compute interpolation coefficients
                      call PDAFomi_get_interp_coeff_lin(4, 2, gcoords, ocoord_p(:, cnt), &
                           thisobs%icoeff_p(:,cnt))

                   end if wetpoint
                endif haveobs
             enddo
          enddo

       else obsmodeB

          ! *** Super-Obbing ***

          ! Allocate process-local index array
          allocate(thisobs%id_obs_p(1, dim_obs_p))

          ! Loop through model grid
          cnt = 0
          do j = idm_start, idm_end
             do i = idm_w, idm_e

                ! Model grid point coordinates
                latM = lat_model(j)
                lonM = lon_model(i)
             
                ! Compute observation grid point indices
                ido_w = ceiling(abs(lonM - lon_obs(1)) / dlonO)
                ido_e = ceiling(abs(lonM - lon_obs(1)) / dlonO)+1
                ido_n = ceiling(abs(latM - lat_obs(1)) / dlatO)+1
                ido_s = ceiling(abs(latM - lat_obs(1)) / dlatO)
             
                ! For wet grid points: Check whether there are valid
                ! observations at the observation grid points around the point
                wetpointB: if (idx_nwet(i, j)>0) then

                   obsflag = 0
                   if (obs_from_file(ido_w,ido_s) > missing_value) obsflag = obsflag + 1
                   if (obs_from_file(ido_e,ido_s) > missing_value) obsflag = obsflag + 1
                   if (obs_from_file(ido_w,ido_n) > missing_value) obsflag = obsflag + 1
                   if (obs_from_file(ido_e,ido_n) > missing_value) obsflag = obsflag + 1

                   ! Use obs. if at least one observation exist in direct vicinity
                   obsflg: if (obsflag>0) then

                      cnt = cnt+1

                      ! Set index of grid point 
                      if (use_wet_state==1 .or. use_wet_state==2) then
                         thisobs%id_obs_p(1, cnt) = idx_nwet(i, j) + sfields(id_obs)%off
                      else
                         thisobs%id_obs_p(1, cnt) = i + nlons*(j-1) + sfields(id_obs)%off
                      end if

                      ! Store observation coordinates
                      ocoord_p(1, cnt) = lon_model(i)
                      ocoord_p(2, cnt) = lat_model(j)

                      ! Compute observation value by averaging
                      obs_sum = 0
                      if (obs_from_file(ido_w,ido_s) > missing_value) &
                           obs_sum = obs_sum + obs_from_file(ido_w,ido_s)
                      if (obs_from_file(ido_e,ido_s) > missing_value) &
                           obs_sum = obs_sum + obs_from_file(ido_e,ido_s)
                      if (obs_from_file(ido_w,ido_n) > missing_value) &
                           obs_sum = obs_sum + obs_from_file(ido_w,ido_n)
                      if (obs_from_file(ido_e,ido_n) > missing_value) &
                           obs_sum = obs_sum + obs_from_file(ido_e,ido_n)
                      obs_p(cnt) = real(obs_sum) / real(obsflag) * sst_scale 

                   end if obsflg

                end if wetpointB

             end do
          end do


       end if obsmodeB

    else

       ! for DIM_OBS_P=0

       allocate(obs_p(1))
       allocate(ivar_obs_p(1))
       allocate(ocoord_p(thisobs%ncoord, 1))
       allocate(thisobs%id_obs_p(1, dim_obs_p))

       if (observation_mode==0) allocate(thisobs%icoeff_p(4, 1))

    end if


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ivar_obs_p(:) = 1.0 / (rms_obs*rms_obs)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! This routine is generic for the case that only the observations, 
    ! inverse variances and observation coordinates are gathered

    call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, lradius_sst_cmems, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!   IF (twin_experiment .AND. filtertype/=11) THEN
!      CALL read_syn_obs(file_syntobs_TYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
!   END IF


! ********************
! *** Finishing up ***
! ********************

    if (doassim_now) then
       ! Deallocate all local arrays
       deallocate(obs_p, ocoord_p, ivar_obs_p)
       deallocate(obs_from_file, lon_obs, lat_obs)
    end if

  end subroutine init_dim_obs_sst_cmems



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  subroutine obs_op_sst_cmems(dim_p, dim_obs, state_p, ostate)

    use PDAFomi, &
         only: PDAFomi_obs_op_interp_lin, PDAFomi_obs_op_gridpoint

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< PE-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    real(pwp), intent(in)    :: state_p(dim_p)        !< PE-local model state
    real(pwp), intent(inout) :: ostate(dim_obs)       !< Full observed state



! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    if (thisobs%doassim==1) then
    
       if (observation_mode==0) then

          ! Observation operator for averaging over grid points
          call PDAFomi_obs_op_interp_lin(thisobs, 4, state_p, ostate)

       else

          ! Observation operator for averaging over grid points
          call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

       end if
    end if

  end subroutine obs_op_sst_cmems



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  subroutine init_dim_obs_l_sst_cmems(domain_p, step, dim_obs, dim_obs_l)

    use PDAFomi, &
         only: PDAFomi_init_dim_obs_l
    use nemo_pdaf, &
         only: wet_pts, nwet
    use assimilation_pdaf, &
         only: domain_coords, locweight

    implicit none

! *** Arguments ***
    integer, intent(in)  :: domain_p     !< Index of current local analysis domain
    integer, intent(in)  :: step         !< Current time step
    integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
    integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector

! *** Local variables ***
    real(pwp) :: coords(2)


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Here one has to specify the coordinates of the local analysis domain
    ! (coords_l) and the localization variables, which can be different for
    ! each observation type and can be made dependent on the index DOMAIN_P.
    ! coords_l should be set in the call-back routine init_dim_l.
    ! coords_l is domain_coords in HBM

    if (trim(dist_obs) == 'gp') then
       coords(1) = wet_pts(6, domain_p)
       coords(2) = wet_pts(7, domain_p)
    else
       coords(1) = domain_coords(1)
       coords(2) = domain_coords(2)
    end if

    call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords, &
         locweight, lradius_sst_cmems, sradius_sst_cmems, dim_obs_l)

  end subroutine init_dim_obs_l_sst_cmems



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! The routine is called in the analysis step of the localized
!! EnKF. It has to apply localization to the two matrices
!! HP and HPH of the analysis step for the module-type
!! observation.
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  subroutine localize_covar_sst_cmems(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    use PDAFomi, only: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
    use assimilation_pdaf, &   
         only: locweight

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< PE-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of observation vector
    real(pwp), intent(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    real(pwp), intent(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    real(pwp), intent(in)    :: coords_p(:,:)         !< Coordinates of state vector elements



! *************************************
! *** Apply covariance localization ***
! *************************************

    ! Here one has to specify the three localization variables
    ! which can be different for each observation type.

    call PDAFomi_localize_covar(thisobs, dim_p, locweight, lradius_sst_cmems, sradius_sst_cmems, &
         coords_p, HP_p, HPH)

  end subroutine localize_covar_sst_cmems

end module obs_sst_cmems_pdafomi
