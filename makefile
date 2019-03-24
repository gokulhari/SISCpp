
BLAS =-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm -lirc -lgfortran
#BLAS = -L /opt/OpenBLAS/lib -lopenblas -lpthread -lm -I /opt/OpenBLAS/include -lgfortran

# Detecting OS from Makefile, code taken from https://gist.github.com/sighingnow/deee806603ec9274fd47

OSFLAG :=
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	OSFLAG += -fopenmp

endif
ifeq ($(UNAME_S),Darwin)
	OSFLAG += -Xpreprecessor -lomp
endif
	UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
	OSFLAG += #-D AMD64
endif
ifneq ($(filter %86,$(UNAME_P)),)
OSFLAG +=   #-D IA32
endif
ifneq ($(filter arm%,$(UNAME_P)),)
	OSFLAG += #-D ARM
endif


sis:
	g++ -c includes/sis.hpp -lfftw3 $(BLAS) -x c++-header

Ex_01:
	g++ examples/EquationSolving/Ex_01.cpp -o bin/Ex_01 -I includes -lfftw3 $(BLAS) -g

Ex_02:
	g++ examples/EquationSolving/Ex_02.cpp -o bin/Ex_02 -I includes -lfftw3 $(BLAS) -g
Ex_03:
	g++ examples/EquationSolving/Ex_03.cpp -o bin/Ex_03 -I includes -lfftw3 $(BLAS) -g
Ex_04:
	g++ examples/EquationSolving/Ex_04.cpp -o bin/Ex_04 -I includes -lfftw3 $(BLAS) -g
Ex_05:
	g++ examples/Eigenvalues/Ex_05.cpp -o bin/Ex_05 -I includes -lfftw3 $(BLAS)
Ex_06:
	g++ examples/Eigenvalues/Ex_06.cpp -o bin/Ex_06 -I includes -lfftw3 $(BLAS)
Ex_07:
	g++ examples/Eigenvalues/Ex_07.cpp -o bin/Ex_07 -I includes -lfftw3 $(BLAS)
Ex_08:
	g++ examples/Eigenvalues/Ex_08.cpp -o bin/Ex_08 -I includes -lfftw3 $(BLAS)
Ex_09:
	g++ examples/Eigenvalues/Ex_09.cpp -o bin/Ex_09 -I includes -lfftw3 $(BLAS)
Ex_10:
	g++ examples/FrequencyResponses/Ex_10.cpp -o bin/Ex_10 -I includes -lfftw3 $(BLAS)
Ex_11:
	g++ examples/FrequencyResponses/Ex_11.cpp -o bin/Ex_11 -I includes -lfftw3 $(BLAS)
Ex_12:
	g++ examples/FrequencyResponses/Ex_12.cpp -o bin/Ex_12 -I includes -lfftw3 $(BLAS)
Ex_13:
	g++ examples/FrequencyResponses/Ex_13.cpp -o bin/Ex_13 -I includes -lfftw3 $(BLAS)
Ex_14:
	g++ examples/Viscoelastic/Ex_14.cpp -o bin/Ex_14 -I includes -lfftw3 $(BLAS)
Ex_15:
	g++ examples/Viscoelastic/Ex_15.cpp -o bin/Ex_15 -I includes -lfftw3 $(BLAS)
Ex_16:
	g++ examples/Viscoelastic/Ex_16.cpp -o bin/Ex_16 -I includes -lfftw3 $(BLAS) -g
Ex_17:
	g++ examples/Viscoelastic/Ex_17.cpp -o bin/Ex_17 -I includes -lfftw3 $(BLAS)
Ex_18:
	g++ examples/Ex_18.cpp -o bin/Ex_18 -I includes -lfftw3 $(BLAS)
Ex_19:
	g++ examples/Ex_19.cpp -o bin/Ex_19 -I includes -lfftw3 $(BLAS) -g
Ex_20:
	g++ examples/Ex_20.cpp -o bin/Ex_20 -I includes -lfftw3 $(BLAS)
Ex_21:
	g++ examples/Ex_21.cpp -o bin/Ex_21 -I includes -lfftw3 $(BLAS)
Ex_22:
	g++ examples/Ex_22.cpp -o bin/Ex_22 -I includes -lfftw3 $(BLAS)
Ex_23:
	g++ examples/Ex_23.cpp -o bin/Ex_23 -I includes -lfftw3 $(BLAS)
Ex_24:
	g++ examples/Ex_24.cpp -o bin/Ex_24 -I includes -lfftw3 $(BLAS)
Ex_25:
	g++ examples/Ex_25.cpp -o bin/Ex_25 -I includes -lfftw3 $(BLAS) $(OSFLAG) -g

Ex_26:
	g++ examples/Ex_26.cpp -o bin/Ex_26 -I includes -lfftw3 $(BLAS) -g
Ex_27:
	g++ examples/Ex_27.cpp -o bin/Ex_27 -I includes -lfftw3 $(BLAS)
Ex_28:
	g++ examples/Ex_28.cpp -o bin/Ex_28 -I includes -lfftw3 $(BLAS)
Ex_29:
	g++ examples/Ex_29.cpp -o bin/Ex_29 -I includes -lfftw3 $(BLAS)
Ex_30:
	g++ examples/Ex_30.cpp -o bin/Ex_30 -I includes -lfftw3 $(BLAS)
Ex_31:
	g++ examples/Ex_31.cpp -o bin/Ex_31 -I includes -lfftw3 $(BLAS)
Ex_32:
	g++ examples/Ex_32.cpp -o bin/Ex_32 -I includes -lfftw3 $(BLAS)
Ex_33:
	g++ examples/Ex_33.cpp -o bin/Ex_33 -I includes -lfftw3 $(BLAS)
format:
	clang-format -i includes/*.hpp examples/*.cpp test/*.cpp paper/paper1/src/*.cpp -style='{BasedOnStyle: llvm, Standard: Cpp03}'

try_spqr:
	clang++ test/try_spqr.cpp -o bin/try_spqr -I includes -I /usr/local/include/eigen3

try_eigen:
	clang++ test/try_eigen.cpp -o bin/try_eigen -I includes -I /usr/local/include/eigen3
matexp:
	clang++ test/matexp.cpp -o bin/matexp -I includes -I /usr/local/include/eigen3 -I FEAST/3.0/include -L FEAST/3.0/lib/x64 -lfeast_dense -lfeast -I includes -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm -lirc -lgfortran
main:
	g++ test/main.cpp -lfftw3 -o bin/main -I includes

Chebfun_test:
	g++ test/Chebfun_test.cpp -lfftw3 -o bin/Chebfun_test -I includes

docs:
	doxygen config_doxygen

Bc_test:
	g++ test/Bc_test.cpp -lfftw3 -o bin/Bc_test -I includes

Linop_test:
	g++ test/Linop_test.cpp -lfftw3 -o bin/Linop_test -I includes
Linop_operator_test:
	g++ test/Linop_operator_test.cpp -lfftw3 -o bin/Linop_operator_test -I includes $(BLAS)
Visco_2D_Cylinder:
	g++ test/Visco_2D_Cylinder.cpp -lfftw3 -o bin/Visco_2D_Cylinder -I includes $(BLAS)

Newt_2D_Cylinder:
	g++ test/Newt_2D_Cylinder.cpp -lfftw3 -o bin/Newt_2D_Cylinder -I includes $(BLAS)

Newt_3D_pipe:
	g++ test/Newt_3D_pipe.cpp -lfftw3 -o bin/Newt_3D_pipe -I includes $(BLAS)
Visco_3D_pipe:
	g++ test/Visco_3D_pipe.cpp -lfftw3 -o bin/Visco_3D_pipe -I includes $(BLAS)
Newt_3D_pipe_double:
	g++ test/Newt_3D_pipe_double.cpp -lfftw3 -o bin/Newt_3D_pipe_double -I includes $(BLAS)
ChannelVisco:
	g++ paper/ChannelVisco.cpp -lfftw3  -I includes $(BLAS) -o bin/ChannelVisco
# -g -o bin/ChannelVisco
# -L /opt/OpenBLAS/lib -lopenblas -lpthread -lm -I /opt/OpenBLAS/include -g -o bin/ChannelVisco

Channel_visco_vorticity:
	g++ paper/Channel_visco_vorticity.cpp -lfftw3 -o bin/Channel_visco_vorticity -I includes $(BLAS)

test_adjoint:
	g++ test/test_adjoint.cpp -lfftw3 -o bin/test_adjoint -I includes $(BLAS) -g

2DViscoSvd:
	g++ test/2DViscoSvd.cpp -o bin/2DViscoSvd -I includes $(BLAS) -g
2DViscoSvd2:
	g++ test/2DViscoSvd2.cpp -o bin/2DViscoSvd2 -I includes $(BLAS) -g

test_sis_feast:
	clang++ test/test_sis_feast.cpp -o bin/test_sis_feast -I includes -I /usr/local/include/eigen3 -I FEAST/3.0/include -L FEAST/3.0/lib/x64 -lfeast_dense -lfeast -I includes -lopenblas -L /opt/OpenBLAS/lib -I /opt/OpenBLAS/include -lgfortran -lfftw3 $(BLAS)
#-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm -lirc -lgfortran
try_feast:
	clang++ test/try_feast.cpp -o bin/try_feast -I includes -I /usr/local/include/eigen3 -I FEAST/3.0/include -L FEAST/3.0/lib/x64 -lfeast_dense -lfeast -I includes -lopenblas -L /opt/OpenBLAS/lib -I /opt/OpenBLAS/include -lgfortran -lfftw3 $(BLAS)
zgeev:
	clang++ test/zgeev.cpp -o bin/zgeev -I includes -I /usr/local/include/eigen3 -lfftw3 $(BLAS)

zggev:
	clang++ test/zggev.cpp -o bin/zggev -I includes -I /usr/local/include/eigen3 -llapacke -lopenblas -L /opt/OpenBLAS/lib -I /opt/OpenBLAS/include -lgfortran -lfftw3 $(BLAS) -v

linop_composition_test:
	g++ test/linop_composition_test.cpp -o bin/linop_composition_test -I includes -lfftw3 $(BLAS) -g
all:
	make Ex_01 Ex_02 Ex_03 Ex_04 Ex_05 Ex_06 Ex_07 Ex_08 Ex_09 Ex_10 Ex_11 Ex_12 Ex_13 Ex_14 Ex_15 Ex_16 Ex_17 Ex_18 Ex_19 Ex_20 Ex_21 Ex_22 Ex_23
.PHONY: run_all
run_all:
	./bin/Ex_01 && ./bin/Ex_02 &&  ./bin/Ex_03 && ./bin/Ex_04 && ./bin/Ex_05 && ./bin/Ex_06 && ./bin/Ex_07 && ./bin/Ex_08 && ./bin/Ex_09 && ./bin/Ex_10 && ./bin/Ex_11 && ./bin/Ex_12 && ./bin/Ex_13 && ./bin/Ex_14 ./bin/Ex_15 && ./bin/Ex_16 && ./bin/Ex_17 && ./bin/Ex_18 && ./bin/Ex_19 && ./bin/Ex_20 && ./bin/Ex_21 && ./bin/Ex_22 && ./bin/Ex_23
