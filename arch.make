export FFLAGS=-traceback -debug -O2 -static_intel
export LDFLAGS=-L/home/user/REPOSITORY/spglib/lib -lsymspg
export MPIFC=mpif90
MKL=$(MKLROOT)/lib/em64t/libmkl_lapack95_lp64.a -Wl,--start-group	\
$(MKLROOT)/lib/em64t/libmkl_intel_lp64.a				\
 $(MKLROOT)/lib/em64t/libmkl_sequential.a				\
 $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread -lm
export LAPACK=$(MKL)
export LIBS=$(LAPACK)