# Makefile
# GNU Makefile for Voronoi tessellation and sparse phase-field grain growth
# Questions/comments to kellet@rpi.edu (Trevor Keller)

# includes
incdir = include
algodir = algorithms

# IBM XL compiler
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_PATH)/include
BG_LIB = -L$(BG_PATH)/lib

# compilers/flags
#compiler = g++ -O3 -Wall
# if you run on SCOREC, make sure you module load gcc/4.9.2 first to have a g++ versin high enough to support C++11
#compiler = g++ -std=c++11 -O3 -W 
# if you add -O3 to the compiler, you may get weird errors std::bad_alloc(), I guess the compiler optimization did something bad...
compiler = g++ -std=c++11 -W                                                                                                                                                                                                                  
pcompiler = mpic++ -O3 -std=c++0x -Wall
flags = -I$(incdir) -I$(algodir) -I$(algodir)/topology

# RPI CCI AMOS compilers/flags
#qcompiler = mpic++ -g -qarch=qp -qtune=qp -qflag=w -qstrict -qreport
qcompiler = mpic++ -O5 -qarch=qp -qtune=qp -qflag=w -qstrict -qprefetch=aggressive -qsimd=auto -qhot=fastmath -qinline=level=10
#qflags = $(CFLAGS) $(BG_INC) $(BG_LIB) $(LDFLAGS) $(flags)

# ONLY uncomment the following if <module load experimental/zlib> FAILS.
qflags = $(BG_INC) $(BG_LIB) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib

# dependencies
core = $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp \
       $(incdir)/MMSP.sparse.hpp

# the program
graingrowth.out: main.cpp graingrowth.cpp tessellate.hpp $(core)
	$(compiler) -DPHASEFIELD $(flags) $< -o graingrowth.out -lz -pthread

mc: main.cpp graingrowth_MC.cpp tessellate.hpp $(core)
	$(compiler) $(flags) $< -o graingrowth.out -lz -pthread

parallel: main.cpp graingrowth.cpp tessellate.hpp $(core)
	$(pcompiler) -DPHASEFIELD $(flags) -include mpi.h $< -o parallel_GG.out -lz

parallelmc: main.cpp graingrowth_MC.cpp tessellate.hpp $(core)
	$(pcompiler) $(flags) -include mpi.h $< -o parallel_MC.out -lz

bgq: main.cpp graingrowth.cpp tessellate.hpp $(core)
	$(qcompiler) $(qflags) -DBGQ -DSILENT -DPHASEFIELD -DRAW $< -o q_GG.out

bgqmc: main.cpp graingrowth_MC.cpp tessellate.hpp $(core)
	$(qcompiler) $(qflags) -DBGQ -DSILENT -DRAW $< -o q_MC.out

wrongendian: wrongendian.cpp
	$(compiler) $< -o $@.out -lz -pthread

mmsp2vtk: TKmmsp2vti.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2vtkFusionline: mmsp2vtkFusionline.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2vtkHaz: mmsp2vtkHaz.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2xyz: mmsp2xyz.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

CountGrain: CountGrain.cpp $(core) 
	$(compiler) $(flags) $< -o $@ -lz

mmsp2xyzStat: mmsp2xyzStat.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2vtkRecolor: mmsp2vtkRecolor.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2vtikeeporder: mmsp2vti_keeporder.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2vti200orient: mmsp2vti_200orient.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

clean:
	rm -rf graingrowth.out parallel_GG.out parallel_MC.out q_GG.out q_MC.out wrongendian.out
