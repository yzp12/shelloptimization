EIGEN_PATH = /u/yzp12/.local/include/eigen3
LIBIGL_PATH = /scratch/cluster/yzp12/projects/domain_regularized_learner/libigl/include
MPIR_INCLUDE = /scratch/cluster/yzp12/mpir-3.0.0/include
MPFR_INCLUDE = /scratch/cluster/yzp12/mpfr-3.1.6_newton/src
MPFRCXX_INCLUDE = mpfrc++-3.6.2
SUITESPARSE_INCLUDE = /scratch/cluster/yzp12/SuiteSparse/include
ALGLIB = alglib/src
#GMP_INCLUDE = /stage/public/linux/include
#GMP_INCLUDE = /stage/public/share/gmp/gmp-4.1.4/
GMP_INCLUDE = /scratch/cluster/yzp12/gmp-6.1.2

MPIR_LIB = /scratch/cluster/yzp12/mpir-3.0.0/.libs
MPFR_LIB = /scratch/cluster/yzp12/mpfr-3.1.6_newton/src/.libs
CHOLMOD_LIB = /scratch/cluster/yzp12/SuiteSparse/CHOLMOD/Lib
GMP_LIB = /scratch/cluster/yzp12/gmp-6.1.2/.libs

CC = g++

CFLAGS = -std=c++11 -Wall -g -O3 -fno-tree-vectorize -DNDEBUG -D MPREAL_HAVE_EXPLICIT_CONVERTERS -pthread

INCLUDES = -I${GMP_INCLUDE} -I${EIGEN_PATH} -I${LIBIGL_PATH} -I${MPIR_INCLUDE} -I${MPFR_INCLUDE} -I${MPFRCXX_INCLUDE} -I${SUITESPARSE_INCLUDE} -I${ALGLIB}

LFLAGS = -L${MPIR_LIB} -L${MPFR_LIB} -L${CHOLMOD_LIB} -Lcollisiondetection/bin

LIBS = 
SPECLIBS = ${CHOLMOD_LIB}/libcholmod.a ${MPIR_LIB}/libmpir.a collisiondetection/bin/libcollisions.a ${MPFR_LIB}/libmpfr.a ${GMP_LIB}/libgmp.a
SHELLS2 = Shells2

MAINSRCS = ${SHELLS2}/CTCD.cpp ${SHELLS2}/Geometry.cpp ${SHELLS2}/main.cpp ${SHELLS2}/RenderingMesh.cpp ${SHELLS2}/ShellEnergy.cpp ${SHELLS2}/SimulationMesh.cpp
MAINOBJ = $(MAINSRCS:.cpp=.o)
ALGLIBSRCS = ${ALGLIB}/alglibinternal.cpp ${ALGLIB}/alglibmisc.cpp ${ALGLIB}/ap.cpp ${ALGLIB}/dataanalysis.cpp ${ALGLIB}/diffequations.cpp ${ALGLIB}/fasttransforms.cpp ${ALGLIB}/integration.cpp ${ALGLIB}/interpolation.cpp ${ALGLIB}/linalg.cpp ${ALGLIB}/optimization.cpp ${ALGLIB}/solvers.cpp ${ALGLIB}/specialfunctions.cpp ${ALGLIB}/statistics.cpp
ALGLIBOBJ = $(ALGLIBSRCS:.cpp=.o)

MAIN = shelloptimization

.PHONY: clean
all: $(MAIN)

$(MAIN) : $(MAINOBJ) $(ALGLIBOBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) $(LIBS) $(MAINOBJ) $(ALGLIBOBJ) -o $(MAIN) $(SPECLIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	rm shelloptimization
	rm Shells2/*.o
	rm $(ALGLIB)/*.o
