# The executable name
# ===================
EXEC = CAPOW

# Choose the machine you are running on.
# ==========================================================================================================
#MACHINE = SCIAMA2
#MACHINE = RAIJIN
#MACHINE = DOGMATIX
#MACHINE = GETAFIX
#MACHINE = TINAROO
MACHINE = LAPTOP
#MACHINE = CORI
#MACHINE = OZSTAR
#MACHINE = PERLMUTTER
#MACHINE = NT

# Options for optimization
# ========================
OPTIMIZE  = -O3

# Setup libraries
# Here is where you'll need to add the correct filepaths for the libraries
# ========================================================================
ifeq ($(MACHINE),SCIAMA2)
  CC = mpicc
ifdef SINGLE_PRECISION
  FFTW_INCL = -I/opt/gridware/pkg/libs/fftw3_float/3.3.3/gcc-4.4.7+openmpi-1.8.1/include/
  FFTW_LIBS = -L/opt/gridware/pkg/libs/fftw3_float/3.3.3/gcc-4.4.7+openmpi-1.8.1/lib/ -lfftw3f_mpi -lfftw3f
else
  FFTW_INCL = -I/opt/gridware/pkg/libs/fftw3_double/3.3.3/gcc-4.4.7+openmpi-1.8.1/include/
  FFTW_LIBS = -L/opt/gridware/pkg/libs/fftw3_double/3.3.3/gcc-4.4.7+openmpi-1.8.1/lib/ -lfftw3_mpi -lfftw3
endif
  GSL_INCL  = -I/opt/apps/libs/gsl/1.16/gcc-4.4.7/include/
  GSL_LIBS  = -L/opt/apps/libs/gsl/1.16/gcc-4.4.7/lib/  -lgsl -lgslcblas
  MPI_INCL  = -I/opt/gridware/pkg/mpi/openmpi/1.8.1/gcc-4.4.7/include
  MPI_LIBS  = -L/opt/gridware/pkg/mpi/openmpi/1.8.1/gcc-4.4.7/lib/ -lmpi
endif

ifeq ($(MACHINE),RAIJIN)
  CC = mpicc
  FFTW_INCL = -I/apps/fftw3/3.3.5/include/
  FFTW_LIBS = -L/apps/fftw3/3.3.5/lib/ -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/apps/gsl/2.3/include/ 
  GSL_LIBS  = -L/apps/gsl/2.3/lib/ -lgsl -lgslcblas
  MPI_INCL  = -I/apps/openmpi/1.6.3/include/ 
  MPI_LIBS  = -L/apps/openmpi/1.6.3/lib/  -lmpi
endif

ifeq ($(MACHINE),GETAFIX)
  CC = mpicc
  FFTW_INCL = -I/opt/fftw/3.3.6/gnu/openmpi2_eth/include/
  FFTW_LIBS = -L/opt/fftw/3.3.6/gnu/openmpi2_eth/lib/ -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/opt/gsl/2.1/gnu/include/
  GSL_LIBS  = -L/opt/gsl/2.1/gnu/lib/ -lgsl -lgslcblas
  MPI_INCL  = -I/opt/openmpi2/gnu/eth/include
  MPI_LIBS  = -L/opt/openmpi2/gnu/eth/lib  -lmpi
endif

ifeq ($(MACHINE),TINAROO)
  CC = mpicc
  FFTW_INCL = -I/opt/fftw/3.3.6/gnu/openmpi_ib/include/
  FFTW_LIBS = -L/opt/fftw/3.3.6/gnu/openmpi_ib/lib/ -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/opt/gsl/2.1/gnu/include/
  GSL_LIBS  = -L/opt/gsl/2.1/gnu/lib/ -lgsl -lgslcblas
  MPI_INCL  = -I/opt/openmpi2/gnu/ib/include
  MPI_LIBS  = -L/opt/openmpi2/gnu/ib/lib  -lmpi
endif

ifeq ($(MACHINE),LAPTOP)
  CC = /opt/homebrew/bin/mpicc
  FFTW_INCL = -I/opt/homebrew/include/
  FFTW_LIBS = -L/opt/homebrew/lib/ -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/opt/homebrew/include/gsl/
  GSL_LIBS  = -L/opt/homebrew/lib/ -lgsl -lgslcblas
  MPI_INCL  = -I/opt/homebrew/include/openmpi/
  MPI_LIBS  = -L/opt/homebrew/lib/openmpi/  -lmpi
endif

ifeq ($(MACHINE),CORI)
  CC = cc
  FFTW_INCL = -I/opt/cray/pe/fftw/3.3.8.2/mic_knl/include/
  FFTW_LIBS = -L/opt/cray/pe/fftw/3.3.8.2/mic_knl/lib/ -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/global/common/sw/cray/cnl7/haswell/gsl/2.5/intel/19.0.3.199/7twqxxq/include/
  GSL_LIBS  = -L/global/common/sw/cray/cnl7/haswell/gsl/2.5/intel/19.0.3.199/7twqxxq/lib/ -lgsl -lgslcblas
  MPI_INCL  = -I/opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/include/
  MPI_LIBS  = -L/opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/lib/  -lmpich
endif

ifeq ($(MACHINE),OZSTAR)
  CC = mpicc
  FFTW_INCL = -I/apps/skylake/software/mpi/gcc/8.2.0/openmpi/4.0.0/fftw/3.3.8/include/
  FFTW_LIBS = -L/apps/skylake/software/mpi/gcc/8.2.0/openmpi/4.0.0/fftw/3.3.8/lib -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/apps/skylake/software/compiler/gcc/8.2.0/gsl/2.5/include/
  GSL_LIBS  = -L/apps/skylake/software/compiler/gcc/8.2.0/gsl/2.5/lib/ -lgsl -lgslcblas
  MPI_INCL = -I/apps/skylake/software/compiler/gcc/8.2.0/openmpi/4.0.0/include/
  MPI_LIBS = -L/apps/skylake/software/compiler/gcc/8.2.0/openmpi/4.0.0/lib -lmpi
endif

ifeq ($(MACHINE),PERLMUTTER)
  CC = cc
  FFTW_INCL = -I/opt/cray/pe/fftw/3.3.10.6/haswell/include/
  FFTW_LIBS = -L/opt/cray/pe/fftw/3.3.10.6/haswell/lib/ -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/usr/include/
  GSL_LIBS  = -L/usr/lib64/ -lgsl -lgslcblas
  MPI_INCL  = -I/opt/cray/pe/mpich/8.1.28/ofi/gnu/12.3/include/
  MPI_LIBS  = -L/opt/cray/pe/mpich/8.1.28/ofi/gnu/12.3/lib/  -lmpich
  FITS_INCL = -I/global/homes/c/chowlett/codes/cfitsio-4.4.0/include/
  FITS_LIBS = -L/global/homes/c/chowlett/codes/cfitsio-4.4.0/lib/ -lcfitsio
endif

ifeq ($(MACHINE),NT)
  CC = mpicc
  FFTW_INCL = -I/apps/modules/software/FFTW.MPI/3.3.10-gompi-2023a/include/
  FFTW_LIBS = -L/apps/modules/software/FFTW.MPI/3.3.10-gompi-2023a/lib -lfftw3_mpi -lfftw3
  GSL_INCL  = -I/apps/modules/software/GSL/2.7-GCC-11.3.0/include/
  GSL_LIBS  = -L/apps/modules/software/GSL/2.7-GCC-11.3.0/lib -lgsl -lgslcblas
  MPI_INCL = -I/apps/modules/software/OpenMPI/4.1.4-GCC-11.3.0/include/
  MPI_LIBS = -L/apps/modules/software/OpenMPI/4.1.4-GCC-11.3.0/lib -lmpi
endif

# Compile the code
# ================
LIBS   =   -lm $(MPI_LIBS) $(FFTW_LIBS) $(GSL_LIBS)

CFLAGS =   $(OPTIMIZE) $(FFTW_INCL) $(GSL_INCL) $(MPI_INCL) $(OPTIONS)

OBJS   = src/main.o src/vars.o src/read_params.o src/read_data.o

INCL   = src/vars.h src/proto.h Makefile

all: $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f src/*.o src/*~ *~ $(EXEC)
