CXX	:= g++
LINK	:= g++
OMP 	:= -fopenmp

MPLAPACKDIR := $(HOME)/opt/MPLAPACK

CPPFLAGS= -Wall -O3 $(OMP) $(MACRO) -std=c++17 -fpermissive -I$(MPLAPACKDIR)/include/mplapack

LFLAGS		:= $(OMP) -llapack -lblas -lm
LFLAGS_BLAS     := $(OMP) -static -llapack -lblas -lgfortran -lquadmath -lm -L/usr/lib/x86_64-linux-gnu/blas -L/usr/lib/x86_64-linux-gnu/lapack
LFLAGS_OPENBLAS := $(OMP) -static -llapack -lblas -lgfortran -lquadmath -lm -L/usr/lib/x86_64-linux-gnu/openblas-openmp

LFLAGS_MPLAPACK1 := $(OMP) -lmplapack__Float64x -lmpblas__Float64x -L$(MPLAPACKDIR)/lib -lopenblas
LFLAGS_MPLAPACK2 := $(OMP) -lmplapack_double -lmpblas_double       -L$(MPLAPACKDIR)/lib -lopenblas
LFLAGS_MPLAPACK3 := $(OMP) -lmplapack_double -lmpblas_double_opt   -L$(MPLAPACKDIR)/lib -lopenblas

EXEFILES := run_lu run_lu_blas run_lu_openblas run_lu_F64 run_lu_double run_lu_double_opt

default: $(EXEFILES)

run_lu: main_lu.o utils.o
	$(LINK) -o $@ $^ $(LFLAGS)

run_lu_blas: main_lu.o utils.o
	$(LINK) -o $@ $^ $(LFLAGS_BLAS)

run_lu_openblas: main_lu.o utils.o
	$(LINK) -o $@ $^ $(LFLAGS_OPENBLAS)

run_lu_F64: main_lu_F64.o utils.o
	$(LINK) -o $@ $^ $(LFLAGS_MPLAPACK1)

run_lu_double: main_lu_double.o utils.o
	$(LINK) -o $@ $^ $(LFLAGS_MPLAPACK2)

run_lu_double_opt: main_lu_double.o utils.o
	$(LINK) -o $@ $^ $(LFLAGS_MPLAPACK3)

clean:;
	rm -rf run_lu $(EXEFILES) *.o

