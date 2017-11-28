EIGEN_PATH = /u/ziwei95/DeepCloth/eigen
LIBIGL_PATH = /u/ziwei95/DeepCloth/libigl/include
MPIR_INCLUDE = /u/ziwei95/mpir-3.0.0/include
MPFR_INCLUDE = /u/ziwei95/mpfr-3.1.6/src
MPFRCXX_INCLUDE = mpfrc++-3.6.2
SUITESPARSE_INCLUDE = /u/ziwei95/SuiteSparse/include
ALGLIB = alglib/src

MPIR_LIB = /u/ziwei95/mpir-3.0.0/.libs
MPFR_LIB = /u/ziwei95/mpfr-3.1.6/src/.libs
CHOLMOD_LIB = /u/ziwei95/SuiteSparse/CHOLMOD/Lib
THREAD_LIB = /usr/lib/x86_64-linux-gnu

CC = g++

CFLAGS = -std=c++11 -Wall -g -O3 -DNDEBUG -D MPREAL_HAVE_EXPLICIT_CONVERTERS -pthread

INCLUDES = -I${EIGEN_PATH} -I${LIBIGL_PATH} -I${MPIR_INCLUDE} -I${MPFR_INCLUDE} -I${MPFRCXX_INCLUDE} -I${SUITESPARSE_INCLUDE} -I${ALGLIB}

LFLAGS = -L${MPIR_LIB} -L${MPFR_LIB} -L${CHOLMOD_LIB} -Lcollisiondetection/bin

LIBS = -lcholmod -lmpir
SPECLIBS = collisiondetection/bin/libcollisions.a ${MPFR_LIB}/libmpfr.a ${THREAD_LIB}/libgmp.a
SHELLS2 = Shells2

MAINSRCS = ${SHELLS2}/CTCD.cpp ${SHELLS2}/Geometry.cpp ${SHELLS2}/main.cpp ${SHELLS2}/RenderingMesh.cpp ${SHELLS2}/ShellEnergy.cpp ${SHELLS2}/SimulationMesh.cpp
MAINOBJ = $(MAINSRCS:.cpp=.o)
ALGLIBSRCS = ${ALGLIB}/alglibinternal.cpp ${ALGLIB}/alglibmisc.cpp ${ALGLIB}/ap.cpp ${ALGLIB}/dataanalysis.cpp ${ALGLIB}/diffequations.cpp ${ALGLIB}/fasttransforms.cpp ${ALGLIB}/integration.cpp ${ALGLIB}/interpolation.cpp ${ALGLIB}/linalg.cpp ${ALGLIB}/optimization.cpp ${ALGLIB}/solvers.cpp ${ALGLIB}/specialfunctions.cpp ${ALGLIB}/statistics.cpp
ALGLIBOBJ = $(ALGLIBSRCS:.cpp=.o)

MAIN = shelloptimization

all: $(MAIN)

$(MAIN) : $(MAINOBJ) $(ALGLIBOBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) $(LIBS) $(MAINOBJ) $(ALGLIBOBJ) -o $(MAIN) $(SPECLIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@
