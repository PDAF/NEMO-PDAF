#!/bin/bash
#SBATCH --job-name=NPDA
#SBATCH -p mpp
#SBATCH --nodes=2
#SBATCH --tasks-per-node=96
#SBATCH --time=0:30:00
#SBATCH --partition=standard96:test
##SBATCH --constraint=turbo_off

# Run script for HLRN-EMMY using Intel compiler version 2022

# To prepare only the run directories without running
# set $prepare=1, $dorun=0, $postproc=0 
# and run interactively in the shell

#set -xv #debugging

module load intel/2022.2
module load openmpi/intel/4.1.4

ulimit -s unlimited

export LD_LIBRARY_PATH=/sw/dataformats/netcdf-parallel/ompi/intel.22/4.9.1/skl/lib:/sw/dataformats/hdf5-parallel/ompi/intel.22/1.12.1/skl/lib:/sw/tools/oneapi/2022.2/compiler/2022.1.0/linux/compiler/lib/intel64_lin:/sw/tools/oneapi/2022.2/mkl/2022.1.0/lib/intel64

# ---------------------------------------------------------------------------------------------------

yy=0001 #Start year
mm=01   #Start month
dd=01   #Start day

tstr_end='00010131' #End date yyyymmdd

np_nemo=46  #number of MPI tasks for each nemo
np_xios=2   #number of MPI tasks for each xios

NENS=4       # Ensemble size

# Whether NEMO should write its outputs (i.e. separate files for each ensemble task)
nemooutput1=0         # (1) Let NEMO write for ensemble member 1
nemooutput_ens=0      # (1) Let NEMO write for ensemble members 2-N
piscesoutput1=0       # (1) Let NEMO-PISCES write for ensemble member 1
piscesoutput_ens=0    # (1) Let NEMO-PISCES write for ensemble members 2-N

# Whether the script prepares the run directories, runs the experiment, does posptprocessing
prepare=1
dorun=1
postproc=1

# Time step size
rn_rdt=5400

# ---------------------------------------------------------------------------------------------------

# Set initial date  $yy$mm$dd==$initial_date we use the restart files from the input directory

initial_date=00010101

# ---------------------------------------------------------------------------------------------------
# Name of experiment and output directory
cn_exp='ORCA2'      # Name of NEMO configuration
EXP='test_DA'       # Name of output directory for storing outputs in post-processing step
# ---------------------------------------------------------------------------------------------------

nemo_exe_dir='/home/hbknerge/SEAMLESS/r4.0.7/cfgs/ORCA2-PDAF_nobio_test/BLD/bin'
xios_exe_dir='/home/hzfblner/SEAMLESS/xios-2.5/bin'
forcing_dir='/scratch/usr/hbknerge/SEAMLESS/run/ORCA2_inputs/inputs'
inputs_nc='/scratch/usr/hbknerge/SEAMLESS/run/ORCA2_inputs/inputs'
setup_store='/scratch/usr/hbknerge/SEAMLESS/run/ORCA2_inputs/config'
initialdir='/scratch/usr/hzfblner/SEAMLESS/restart'
restart_out='output/restarts'
# ---------------------------------------------------------------------------------------------------

#export nemo_exe_dir
#export restart_out

# ---------------------------------------------------------------------------------------------------

tstr_ini="$yy$mm$dd"                          # Date at start of run
tstrm1=$(date -I -d "$tstr_ini - 1 day ")     # Previous day (used to link forcings)
tstrendp1=$(date -I -d "$tstr_end + 1 day ")  # Final day plus one (used to link forcings)
tstr=$tstr_ini                                # Store date at start of run

echo '----------------------------------------'
echo 'NEMO-PDAF data assimilation'
if [ $tstr -eq $initial_date ]; then
    echo "Start from initial date: " $initial_date
else
    echo "Restart at date:         " $tstr
fi
echo "Final date:              " $tstr_end
echo "Experiment:              " $EXP
echo "Results are stored in: " $EXP
echo 'Ensemble size = ' $NENS
echo '----------------------------------------'

# ---------------------------------------------------------------------------------------------------

# Prepare run directories and files
if [ $prepare -eq 1 ]; then

    if [ $tstr -eq $initial_date ]; then

        # Create working directories for different ensemble members
	echo ' '
	echo 'Creating ensemble working directories...'
	for((i=1;i<=$NENS;i++))
	  do
	  ENSstr=`printf %03d $i`
	  if [ ! -d ${ENSstr} ]; then
	      mkdir -p ${ENSstr}

              wdir=`pwd`/${ENSstr}
              mkdir -p $wdir/output/restarts
	      mkdir -p $wdir/output/log
	      mkdir -p $wdir/initialstate
	      mkdir -p $wdir/forcing
          fi

	  echo 'Run directory: ' $wdir
	done

        # Create experiment output directories
        echo  ' '
        echo 'Create experiment output directories...'
        mkdir -p $EXP/cfg
        mkdir -p $EXP/data
        mkdir -p $EXP/DA

	echo ' '

	# Link executables into run directories
	echo 'linking executables...'
	for((i=1;i<=$NENS;i++))
	  do
	  ENSstr=`printf %03d $i`
	  wdir=`pwd`/${ENSstr}
	  if [ ! -f ${ENSstr}/xios_server.exe ]; then
	      ln -s $xios_exe_dir/xios_server.exe $wdir/
	  fi
	  if [ ! -f ${ENSstr}/nemo.exe ]; then
	      ln -s $nemo_exe_dir/nemo.exe $wdir/
	  fi
	done

    fi # if [ $tstr -eq $initial_date ]

# ---------------------------------------------------------------------------------------------------

    echo 'Simulation with nemo executable from ' $nemo_exe_dir
    

# ---------------------------------------------------------------------------------------------------

    # Prepare PDAF namelist
    if [ $tstr -eq $initial_date ]; then
        pdafrestart=".false."
    else
        pdafrestart=".true."
    fi

    cat namelist_cfg.pdaf_template     \
    | sed -e "s:_DIMENS_:$NENS:"     \
    | sed -e "s:_RESTART_:$pdafrestart:"     \
    | sed -e "s:_MONTH_:$mm:"     \
    > namelist_cfg.pdaf

# ---------------------------------------------------------------------------------------------------

    # Linking restart files - optional

    link_restarts=0  # Set to 1 to activate this block
    if [ $link_restarts -eq 1 ]; then
      for((i=1;i<=$NENS;i++)); do
        ENSstr=`printf %03d $i`
        wdir=`pwd`/${ENSstr}

        if [ $i -eq 1 ]; then
	  echo '----------------------------------------'
	  echo 'Prepare restart files'
        fi
        if [ $tstr -eq $initial_date ]; then
          if [ $i -eq 1 ]; then
              echo "Link initial restart files"
          fi
#	  if [ ! -f $wdir/restart_in.nc ]; then
#              ln -s $initialdir/${nc_exp}_restart.nc $wdir/restart_in.nc
#	  fi
        fi
      done
    fi # if [ $link_restarts -eq 1 ]; then

    
# ---------------------------------------------------------------------------------------------------

    echo '----------------------------------------'
    echo 'Preparing forcing'
    echo 'Time period from ' $tstrm1 'until' $tstrendp1

    if [ $tstr -eq $initial_date ]; then

	echo 'Linking files with fixed values ...'
	for((i=1;i<=$NENS;i++))
	  do
	  ENSstr=`printf %03d $i`
	  wdir=`pwd`/${ENSstr}

          # For the ORCA2 test case, we link all input files from a single directory
          cd $wdir
          ln -s $inputs_nc/* .          
          cd -
	done
    fi # if [ $tstr -eq $initial_date ]

    
    # Link time-dependent forcing files - optional
    
    link_time_dep_forcing=0  # Set to 1 to activate this block
    if [ $link_time_dep_forcing -eq 1 ]; then
      tstrhere=$tstrm1
      while [ "$(date -d "$tstrhere" +%Y%m%d)" -le "$(date -d "$tstrendp1" +%Y%m%d)" ]; do
	date_nemo=$(date -d "$tstrhere" +y%Ym%md%d)  # Date format yYYYYmMMdDD
	date_force=$(date -d "$tstrhere" +%Y%m%d)    # Date format YYYYMMDD
	
	tstrhere=$(date -d "$tstrhere" +%Y%m%d)      # Loop date

	echo 'Link time-varying forcing files ... date ' $tstrhere
	for((i=1;i<=$NENS;i++)); do
	  ENSstr=`printf %03d $i`
	  wdir=`pwd`/${ENSstr}

          # Example: Link TS boundary file
	  # This is configures for daily files
#	  if [ -f $wdir/forcing/bdy_ts_$date_nemo.nc ]; then
#	    if [ $i -eq 1 ]; then
#             echo 'Use existing file ' bdy_ts_$date_nemo.nc
#	    fi
#	  else
#	    ln -s $forcing_dir/bdy_ts_$tstrhere.nc  $wdir/forcing/bdy_ts_$date_nemo.nc
#	  fi
	done

	tstrhere=$(date -I -d "$tstrhere + 1 day ")
	
      done
    fi

    echo '----------------------------------------'

    # Count days
    ndays=0
    while [ "$(date -d "$tstr" +%Y%m%d)" -le "$(date -d "$tstr_end" +%Y%m%d)" ]; do
	date_nemo=$(date -d "$tstr" +y%Ym%md%d)  # Date format yYYYYmMMdDD
	tstr=$(date -I -d "$tstr + 1 day ")
	ndays=$(($ndays + 1))
    done
    echo 'Number of days in this run: ' $ndays

    # Compute number of time steps in this run
    nn_itend=$(($ndays * 24 * 3600 / $rn_rdt ))
    nnstep=`printf %08d $nn_itend`
    echo 'run length nn_itend: ' $nn_itend

    echo '----------------------------------------'

# ---------------------------------------------------------------------------------------------------

    # Prepare XML files for ensemble outputs
    echo 'Prepare XML files ...'
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      ENSstr2=`printf %02d $i`

      cat $setup_store/context_nemoXXX.xml_template   \
	  | sed -e "s:_DIMENS3_:${ENSstr}:g"   \
	  > context_nemo${ENSstr}.xml

      cat $setup_store/file_def_nemo-iceXXX.xml_template   \
	  | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
	  > file_def_nemo-ice${ENSstr}.xml

      cat $setup_store/file_def_nemo-oceXXX.xml_template   \
	  | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
	  > file_def_nemo-oce${ENSstr}.xml

      cat $setup_store/file_def_nemo-piscesXXX.xml_template   \
          | sed -e "{s:_DIMENS3_:${ENSstr}:g;s:_DIMENS2_:${ENSstr2}:g}"   \
          > file_def_nemo-pisces${ENSstr}.xml
    done

    cp  $setup_store/iodef.xml_template iodef.xml
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      echo "  <context id=\"nemo_${ENSstr}\" src=\"./context_nemo${ENSstr}.xml\"/> " >> iodef.xml
    done
    echo "</simulation>" >> iodef.xml

    #define NEMO_001 output - ensemble member 1
    if [ $nemooutput1 -eq 1 ]; then
	FILE101flag=.true.
    else
	FILE101flag=.false.
	echo "Deactivate NEMO output for ensemble member 1"
    fi

    #define NEMO_00X output - ALL ENSEMBLE MEMBERS 2-n
    if [ $nemooutput_ens -eq 1 ]; then
	FILE10Xflag=.true.
    else
	FILE10Xflag=.false.
	echo "Deactivate NEMO output for ensemble members 2-NENS"
    fi

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      if [ $i -eq 1 ]; then
	  sed -i "s:_FILE1flag_:${FILE101flag}:g" file_def_nemo-oce001.xml
	  sed -i "s:_FILE1flag_:${FILE101flag}:g" file_def_nemo-ice001.xml
      else
	  sed -i "s:_FILE1flag_:${FILE10Xflag}:g" file_def_nemo-oce${ENSstr}.xml
	  sed -i "s:_FILE1flag_:${FILE10Xflag}:g" file_def_nemo-ice${ENSstr}.xml
      fi
    done

    #define PISCES_001 output - ensemble member 1
    if [ $piscesoutput1 -eq 1 ]; then
	FILE201flag=.true.
    else
	FILE201flag=.false.
    fi

    #define PISCES_00X output - ALL ENSEMBLE MEMBERS 2-n
    if [ $piscesoutput_ens -eq 1 ]; then
	FILE20Xflag=.true.
    else
	FILE20Xflag=.false.
    fi

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      if [ $i -eq 1 ]; then
	  sed -i "s:_FILE2flag_:${FILE201flag}:g" file_def_nemo-pisces001.xml
      else
	  sed -i "s:_FILE2flag_:${FILE20Xflag}:g" file_def_nemo-pisces${ENSstr}.xml
      fi
    done

# ---------------------------------------------------------------------------------------------------

    # Prepare NEMO namelist file and copy configuration files
    echo "Copy and prepare configuration files"

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}

      cp $setup_store/*.xml $wdir/
      cp $setup_store/namelist* $wdir/
      cp `pwd`/namelist_cfg.pdaf $wdir/
      cp `pwd`/*.xml $wdir/

      cat $wdir/namelist_cfg_template                  \
          | sed -e "s:_nn_itend_:$nn_itend:"    \
          | sed -e "s:_nn_date0_:$tstr_ini:"        \
          | sed -e "s:_rn_rdt_:$rn_rdt:"        \
          > $wdir/namelist_cfg

    done
    # Clean up xml prepared xml files
    rm `pwd`/*.xml 


# ---------------------------------------------------------------------------------------------------

    # Create MPMD configuration file for srun

    echo 'Prepare mpmd.conf ...'

    if [ -f mpmd.conf ]; then
	rm mpmd.conf
    fi

    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      #echo 'preparing mpmd.conf...'${ENSstr}

      echo '#!/bin/sh' > nemo${ENSstr}.sh
      echo 'cd '${ENSstr} >> nemo${ENSstr}.sh
      echo './nemo.exe' >> nemo${ENSstr}.sh
      chmod +x nemo${ENSstr}.sh

      echo '#!/bin/sh' > xios${ENSstr}.sh
      echo 'cd '${ENSstr} >> xios${ENSstr}.sh
      echo './xios_server.exe' >> xios${ENSstr}.sh
      chmod +x xios${ENSstr}.sh

      echo $(((i-1)*($np_nemo+$np_xios)))'-'$(((i-1)*($np_nemo+$np_xios)+$np_nemo-1))' ./nemo'${ENSstr}.sh >> mpmd.conf
      echo $(((i-1)*($np_nemo+$np_xios)+$np_nemo))'-'$(((i)*($np_nemo+$np_xios)-1))' ./xios001.sh' >> mpmd.conf
    done
    
    cat mpmd.conf
    
    echo 'JOB PREPARATIONS COMPLETED'
    echo '----------------------------------------'

fi  # if $prepare==1


# ---------------------------------------------------------------------------------------------------

# Execute the run
if [ $dorun -eq 1 ]; then
    echo 'RUN ASSIMILATION'

        srun -l --cpu_bind=cores --multi-prog mpmd.conf

    # error check in ocean.output
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}
  
      if grep -q 'E R R O R' $wdir/ocean.output
        then
          { echo 'ERROR in ocean.output' ;  }
        exit 1
      fi
    done

fi # if $dorun==1


# ---------------------------------------------------------------------------------------------------
# POSTPROCESSING

if [ $postproc -eq 1 ]; then

    echo '----------------------------------------'
    echo "POSTPROCESSING"
    date

# ---------------------------------------------------------------------------------------------------
# Move output files

    echo "Move NEMO and DA output files from directory 001"
    wdir=`pwd`

    mv $wdir/001/state_*.nc $EXP/DA/
    mv $wdir/001/variance_*.nc $EXP/DA/
    mv $wdir/001/${cn_exp}_*_???.nc $EXP/data/
    cp $wdir/001/namelist_cfg.pdaf $EXP/cfg/namelist_cfg.pdaf_$tstr_ini-$tstr



# ---------------------------------------------------------------------------------------------------
#Remove old files

    echo "clean up old files"
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}

      # Her you might like to e.g. remove forcing files
#      rm $wdir/forcing/bdy_ts_*.nc
    done


# ---------------------------------------------------------------------------------------------------
#save restarts

    echo "Save restart files"

    np=$(( $np_nemo - 1 ))
    for((i=1;i<=$NENS;i++))
      do
      ENSstr=`printf %03d $i`
      wdir=`pwd`/${ENSstr}

      # Save log file
      cp $wdir/ocean.output $wdir/output/log/ocean.output_${tstr_ini}-${tstr}

      # Rename restart files to include date
      for n in `seq -f "%04g" 0 $np`;do
	  mv $wdir/${cn_exp}'_'$nnstep'_restart_'$n'.nc'     $wdir/$restart_out/'/restart_'$n'_'$date_nemo'.nc'
	  mv $wdir/${cn_exp}'_'$nnstep'_restart_ice_'$n'.nc' $wdir/$restart_out/'/restart_ice_'$n'_'$date_nemo'.nc'
      done
    done

    echo "END POSTPROCESSING"
    date

fi # if postproc=1

# ---------------------------------------------------------------------------------------------------

echo 'JOB DONE'
exit
