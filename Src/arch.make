# arch.make

# Compiler
export CPU_COMPILER = mpif90
export GPU_COMPILER = nvfortran

# CPU and GPU FFLAGS
CPU_FFLAGS = -qopenmp -traceback -O2 -fpp -DCPU_VERSION  #-static_intel   -debug 
GPU_FFLAGS = -acc -Minfo=accel  -Mpreprocess -g -cuda -gpu=ptxinfo -mp -O2 -DGPU_VERSION -DGPU_ALL_MODE_PARALLELIZATION
# -DGPU_ALL_MODE_PARALLELIZATION # use this tag to use all-mode parallelization
# -DGPU_MODE_BY_MODE_PARALLELIZATION # change to this tag to use mode-by-mode parallelization 

# Linking
# Example paths - replace with your actual library paths:
# LDFLAGS = -L/path/to/your/lib -lsymspg
LDFLAGS = -lsymspg
LDFLAGS += -latomic
export LDFLAGS

MKL = $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
      $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
      $(MKLROOT)/lib/intel64/libmkl_sequential.a \
      $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

export LAPACK = $(MKL)
export LIBS = $(LAPACK)

