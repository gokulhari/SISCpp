
BLAS =-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm -lirc -lgfortran
#BLAS = -L /opt/OpenBLAS/lib -lopenblas -lpthread -lm -I /opt/OpenBLAS/include -lgfortran

# Detecting OS from Makefile, code taken from https://gist.github.com/sighingnow/deee806603ec9274fd47
PYTHON = -I/Users/gokul/anaconda3/include/python3.7m -L/Users/gokul/anaconda3/lib -lpython3.7
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
2DViscoSvd2:
	g++ test/2DViscoSvd2.cpp -o bin/2DViscoSvd2 -I includes $(BLAS) -g
2DViscoSvd_LikeMJ:
	g++ test/2DViscoSvd_LikeMJ.cpp -o bin/2DViscoSvd_LikeMJ -I includes $(BLAS) -g
2DViscoSvdOmegaVar:
	g++ test/2DViscoSvdOmegaVar.cpp -o bin/2DViscoSvdOmegaVar -I includes $(BLAS) -g

format:
	clang-format -i includes/*.hpp examples/*.cpp test/*.cpp paper/paper1/src/*.cpp -style='{BasedOnStyle: llvm, Standard: Cpp03}'

docs:
	doxygen config_doxygen

all:
	make Ex_01 Ex_02 Ex_03 Ex_04 Ex_05 Ex_06 Ex_07 Ex_08 Ex_09 Ex_10 Ex_11 Ex_12 Ex_13 Ex_14 Ex_15 Ex_16 Ex_17
.PHONY: run_all
run_all:
	./bin/Ex_01 && ./bin/Ex_02 &&  ./bin/Ex_03 && ./bin/Ex_04 && ./bin/Ex_05 && ./bin/Ex_06 && ./bin/Ex_07 && ./bin/Ex_08 && ./bin/Ex_09 && ./bin/Ex_10 && ./bin/Ex_11 && ./bin/Ex_12 && ./bin/Ex_13 && ./bin/Ex_14 ./bin/Ex_15 && ./bin/Ex_16 && ./bin/Ex_17
