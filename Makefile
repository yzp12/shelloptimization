EIGEN_PATH = /u/yzp12/.local/include/eigen3
#EIGEN_PATH = /u/evouga/eigen
LIBIGL_PATH = /scratch/cluster/yzp12/projects/domain_regularized_learner/libigl/include
#LIBIGL_PATH = /u/evouga/tools/libigl/include
SUITESPARSE_INCLUDE = /scratch/cluster/yzp12/SuiteSparse/include
#SUITESPARSE_INCLUDE = /u/evouga/tools/SuiteSparse/include
ALGLIB = alglib/src

CHOLMOD_LIB = /scratch/cluster/yzp12/SuiteSparse/CHOLMOD/Lib
#CHOLMOD_LIB = /u/evouga/tools/SuiteSparse/lib

CC = g++

CFLAGS = -std=c++11 -Wall -g -O3 -fno-tree-vectorize -DNDEBUG -D MPREAL_HAVE_EXPLICIT_CONVERTERS -pthread

INCLUDES = -I${EIGEN_PATH} -I${LIBIGL_PATH} -I${SUITESPARSE_INCLUDE} -I${ALGLIB}

LFLAGS = -L${CHOLMOD_LIB} -Lcollisiondetection/bin

LIBS = 
SPECLIBS = -L${CHOLMOD_LIB} -lcholmod collisiondetection/bin/libcollisions.a
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
