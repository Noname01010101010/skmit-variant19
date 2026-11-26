CXX     ?= g++
MPICXX  ?= mpicxx
OMPFLAG ?= -fopenmp
CXXFLAGS = -O1 -std=c++11

all: seq omp mpi hybrid

seq: seq.cpp common.h
	$(CXX) $(CXXFLAGS) seq.cpp -o seq
omp: omp.cpp common.h
	$(CXX) $(CXXFLAGS) $(OMPFLAG) omp.cpp -o omp
mpi: mpi.cpp common.h partition.h
	$(MPICXX) $(CXXFLAGS) mpi.cpp -o mpi
hybrid: hybrid.cpp common.h partition.h
	$(MPICXX) $(CXXFLAGS) $(OMPFLAG) hybrid.cpp -o hybrid

test:
	bash tests/run_tests.sh

clean:
	rm -rf seq omp mpi hybrid tests/build

.PHONY: all test clean
