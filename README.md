# NEMO-PDAF

NEMO-PDAF coupling routines. For the contributors to the code see MY_SRC_pdaf/initialize_pdaf.F90

For documentation on how to use the code, see the Wiki on https://github.com/PDAF/NEMO-PDAF.

Directories:

**MY_SRC_pdaf/** contains the pdaf user routines

**MY_SRC_nemo4.0/** contains the modified NEMO routines for NEMO 4.0 (e.g. NEMO-NORDIC)

**MY_SRC_NEMO4.0.7/** contains the modified NEMO routines for NEMO 4.0.7. This should also be usable for other 4.0.x versions.

**run/** contains 
- example `arch` file for compiling NEMO with PDAF (arch-linux_mpifort_PDAF.fcm)
- namelist file (namelist_cfg.pdaf_template)
- example run script for SLURM (run_NEMO-PDAF_N4.sl)

**run/config/** contains template files for configuring XIOS ensemble file output. These files are used by the example run script. 

**TOOLS/generate_covar/** contains the code for the program generating a covariance matrix from model snapshots using a singular value decomposition

**examples/ORCA2_PISCES/** contains the files to run the NEMO test case ORCA2_ICE_PISCES
