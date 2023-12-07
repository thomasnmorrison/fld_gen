#-*-Makefile-*-

FC = gfortran
FOPT = -o2 -fopenmp
FFLAGS = -fdefault-real-8 -fdefault-double-8 -cpp -ffree-line-length-none -Wall

FFTW_INC = -I /usr/include
FFTW_LIB = -L /usr/lib/x86_64-linux-gnu
BLAS_LIB = -L /usr/lib/x86_64-linux-gnu/blas
LAPACK_LIB = -L /usr/lib/x86_64-linux-gnu/lapack
LIBS = -lfftw3_omp -lfftw3 -llapack -lblas -lm

CORR = butcher_table.o gl_integrator.o newton_root.o util.o potential.o eom_bg_cosmic.o eom_pert_cosmic.o corr_cosmic.o corr_mod.o
OBJS = fftw_mod.o params.o vars.o grv_mod.o lin_tran.o nlin_tran.o io_mod.o spec_init.o $(CORR)

executable: %: $(OBJS) fld_gen.o
	$(FC) $(FFLAGS) $(FOPT) $(FFTW_INC) -o genmock fld_gen.f90 $(OBJS) $(FFTW_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) $(FOPT) $(FFTW_INC) -c $< -o $@ $(FFTW_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f genmock
