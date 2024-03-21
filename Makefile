# Export MKLROOT in your session to Intel MKL's installation folder. 
# For example:
# export MKLROOT=${HOME}/usr/intel/mkl
MPIF90 = mpifort
COMPILER_ = icx
FCFLAGS = -diag-disable=10448 -c -O3 -fpp -heap-arrays -mcmodel=medium -fp-model consistent -fpconstant -I${USRLOCAL}/include -I${MKLROOT}/include/mkl/intel64/lp64 -I${MKLROOT}/include -qmkl=sequential
LDFLAGS = -L${USRLOCAL}/lib -lfftw3 -L${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64/ -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# Modules
MODULES = numbers.o\
		  openmpi.o\
		  io.o\
		  parameters.o\
		  fftw.o\
		  diffops.o\
		  symmops.o\
		  test_symmops.o\
		  fieldio.o\
		  vfield.o\
		  rhs.o\
		  timestep.o\
		  stats.o\
		  lyap.o\
		  symred.o\
		  projector.o\
		  run.o\
		  solver.o

# Objects

UTILS 	= utilities/newton.o\
          utilities/eigen.o

# -------------------------------------------------------

all:  $(MODULES) main.f90
	  $(MPIF90) $(FCFLAGS) main.f90
	  $(MPIF90)  -diag-disable=10448 -o dns.x main.o $(MODULES) $(LDFLAGS)
	  $(COMPILER_) -diag-disable=10448 -fPIC -c fakeintel.c
	  $(COMPILER_) -diag-disable=10448 -shared -o libfakeintel.so fakeintel.o

# -------------------------------------------------------

# utils: $(MODULES) $(UTILS)
# 	   $(MPIF90) $(MODULES) newton.o -o newton.x $(LDFLAGS)
# 	   $(MPIF90) $(MODULES) eigen.o -o eigen.x $(LDFLAGS)
	   
# -------------------------------------------------------

# compile

# $(OBJ): $(MODULES) 

%.o: %.f
	$(MPIF90) $(FCFLAGS) $<

m_parameters.o: m_parameters.f90
				bash version.sh && $(MPIF90) $(FCFLAGS) $<

%.o: %.f90
	 $(MPIF90) $(FCFLAGS) $<

clean:
	rm *.o *.mod *.x
