# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################

# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory BASEDIR/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package
BASEDIR = ../../PDAF

# CPP definitions. Define 'key_top' if the state vector contains also biogeochemistry
DEFS = -DPDAF_OFFLINE  #-Dkey_top

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h

# Name of executable
EXE = PDAF_eofcovar

# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =	timer.o \
		mod_memcount_pdaf.o \
		mod_kind_pdaf.o

# PDAF user routines identical to online analysis code
OBJ_PDAF_USER =	parallel_pdaf.o \
		nemo_pdaf.o \
		statevector_pdaf.o \
		transforms_pdaf.o \
		io_pdaf.o

# Routines particular for the offline code
OBJ_OFFLINE =	assimilation_pdaf_eofcovar.o \
		initialize_offline.o \
		main_eofcovar.o \
		eofcovar.o

######################################################

$(EXE) : libpdaf-d.a \
	$(MODULES) $(OBJ_PDAF_USER) $(OBJ_OFFLINE) $(OBJ_MPI)
	$(LD)  $(OPT_LNK) $(NC_LIB) -o $@  \
	$(MODULES) $(OBJ_PDAF_USER) $(OBJ_OFFLINE) $(OBJ_MPI) \
	-L$(BASEDIR)/lib -lpdaf-d \
	$(LINK_LIBS)
	@echo "++++++ Done ++++++"

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR)/src; make;

######################################################

.F90.o :
	$(FC) $(OPT) $(DEFS) $(CPP_DEFS) $(MPI_INC) $(NC_INC) -I$(BASEDIR)/include -c \
	$*.F90

.f90.o :
	$(FC) $(OPT) $(DEFS) $(CPP_DEFS) $(MPI_INC) $(NC_INC) -c \
	$*.f90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

cleanall : cleanpdaf clean

clean :
	rm -f *.o *.mod *.a $(EXE)

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR)/src; make clean


# Some necessary dependencies to ensure that the program compiles
eofcovar.o: assimilation_pdaf_eofcovar.o
io_pdaf.o: assimilation_pdaf_eofcovar.o
