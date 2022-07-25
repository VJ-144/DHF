CXX = g++
#CXX = clang++

# Compiler flags:
#  -DNO_ODE compiles without boost/ode package used for flow equation solver
#  -DNO_HDF5 compiles without boost/ode package used for flow equation solver
#  -DOPENBLAS_NOUSEOMP=1 removes parallel blocks which take threads away from OPENBLAS
#                        to be used if OpenBlas was compiled without the USE_OMP flag
.PHONY: version
ALL = OrbitsTest.exe HFTest.exe
INSTDIR = $(HOME)
INCLUDE   = -I./armadillo
FLAGS     = -O3 -march=native -std=c++11 -fopenmp
#LIBS = -lopenblas -lgsl -lz
LIBS = -lblas -lgsl -lz
PYTHON=off
#DEBUG=on

# I assume we're running on linux
OS = LINUX
# But in case we're crazy enough to run on MacOS, might as well check...
ifneq (,$(findstring arwin,$(shell uname)))
  OS = MACOS
endif

ifeq ($(DEBUG),on)
 FLAGS     = -march=native -std=c++11 -fopenmp -g
endif

WARNFLAGS = -Wall -Wno-comment -Wno-deprecated-declarations -Wno-missing-braces
FLAGS    += $(WARNFLAGS)

SOFLAGS = $(FLAGS)

THEHOST = $(shell if [ `hostname|grep jrl` ]; then echo jureca; elif [ `hostname|grep cougar` ]; then echo cougar; elif [ `hostname|grep cronos` ]; then echo cronos; elif [ `hostname|grep oak` ]; then echo oak; elif [ `hostname|grep cedar` ]; then echo cedar; elif [ `hostname|grep mox` ]; then echo mox; elif [ `hostname|grep bebop` ]; then echo bebop; else echo other; fi)

ifeq ($(OS),MACOS)
  FLAGS     = -Xpreprocessor -fopenmp -O3  -std=c++11 -I/usr/local/include -L/usr/local/lib
  LIBS += -lomp
  PYTHONFLAGS =  $(shell python-config --cflags | sed -e 's/-arch i386//')
endif


ifeq ($(THEHOST),other)  # default options. assumes boost and python are set up nicely.
 LIBS += -llapack
# FLAGS += -DOPENBLAS_NOUSEOMP=1
# ifneq ($(PYTHON),off)
#  ALL += pyIMSRG.so
# endif
endif


all: $(ALL)
#	@echo Building with build version $(BUILDVERSION)

#OBJ = ModelSpace.o TwoBodyME.o ThreeBodyME.o Operator.o  ReadWrite.o\
#      HartreeFock.o imsrg_util.o Generator.o IMSRGSolver.o AngMom.o\
#      IMSRGProfiler.o Commutator.o HFMBPT.o  RPA.o\
#      M0nu.o DarkMatterNREFT.o Jacobi3BME.o UnitTest.o \
#      TwoBodyChannel.o ThreeBodyChannel.o \
#      ThreeBodyStorage.o ThreeBodyStorage_pn.o ThreeBodyStorage_iso.o \
#      ThreeBodyStorage_no2b.o ThreeBodyStorage_mono.o ThreeLegME.o  \
#      boost_src/gzip.o boost_src/zlib.o

boost_src/%.o: boost_src/%.cpp
	$(CXX) -c $^ -o $@ $(INCLUDE) $(FLAGS)

%.o: %.cc %.hh
	$(CXX) -c $*.cc -o $@ $(INCLUDE) $(FLAGS)

%.o: %.cc
	$(CXX) -c $*.cc -o $@ $(INCLUDE) $(FLAGS)

OrbitsTest.exe: OrbitsTest.o Orbits.o TwoBodySpace.o ModelSpace.o Operator.o TwoBodyOperator.o 
	$(CXX) -o $@ $^ $(INCLUDE) $(FLAGS) $(LIBS)
HFTest.exe: HFTest.o HFMBPT.o HartreeFock.o Operator.o TwoBodyOperator.o ModelSpace.o TwoBodySpace.o Orbits.o 
	$(CXX) -o $@ $^ $(INCLUDE) $(FLAGS) $(LIBS)

clean:
	rm -f *.o *.so *.exe boost_src/*.o
