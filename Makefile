#PS_PATH =../fdps/fdps-devel/src/
PS_PATH =../fdps/fdps-test/src/

INC = -I$(PS_PATH)

CC = time g++
#CC = time mpicxx
CFLAGS = -O3 --std=gnu++1z  # 
#CFLAGS += -Wall
CFLAGS += -ffast-math    -ftree-vectorize  -fopt-info-vec-optimized=vector.txt -march=native
CFLAGS +=  -march=native#  -pg
#CFLAGS +=  -mavx2
CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -DPARTICLE_SIMULATOR_USE_SAMPLE_SORT
#CFLAGS += -fsanitize=address -fsanitize=leak -static-libasan
MPICC =  time mpicxx
MPICFLAGS = $(CFLAGS) -DPARTICLE_SIMULATOR_MPI_PARALLEL
SRCS = LICENSE  Readme.md nbody-with-center.cpp user-defined.hpp\
      Makefile ringin	 ring.rb	  samplein \
      Makefile.a64fx Makefile.a64fxclang
EXPORTDIR = ../nbody-with-center-export
use_phantom_grape_x86 = no
#use_gpu_cuda = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

all:nbody-with-center

ifeq ($(use_phantom_grape_x86),yes)
PG_ROOT = $(PS_PATH)/phantom_grape_x86/G5/newton/libpg5
INC += -I$(PG_ROOT)
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   = -L$(PG_ROOT) -lpg5
PG_BUILD = cd $(PG_ROOT) && $(MAKE) distclean libpg5.a
PG_CLEAN = cd $(PG_ROOT) && $(MAKE) distclean
else
PG_BUILD =
PG_CLEAN = 
endif

ifeq ($(use_gpu_cuda),yes)
CUDA_HOME = /usr/local/cuda
#CUDA_HOME = /gwfefs/opt/x86_64/cuda/7.5
NVCC = time $(CUDA_HOME)/bin/nvcc -Xcompiler="-O3"
INC  += -I$(CUDA_HOME)/samples/common/inc/
CFLAGS  += -DENABLE_GPU_CUDA
CLIBS = -L$(CUDA_HOME)/lib64 -lcudart -lgomp
force_gpu_cuda.o:force_gpu_cuda.cu
	$(NVCC) $(INC) -c -o $@ $<
OBJS = force_gpu_cuda.o
endif

nbody-with-center:nbody-with-center.cpp user-defined.hpp  $(OBJS)
	$(PG_BUILD)
	$(CC) $(INC) $(CFLAGS) -pg -o $@ nbody-with-center.cpp $(CLIBS)
nbody-with-center-mpi:nbody-with-center.cpp user-defined.hpp  $(OBJS)
	$(PG_BUILD)
	$(MPICC) $(INC) $(MPICFLAGS) -o $@ nbody-with-center.cpp $(CLIBS)
nbody-with-center-quad-mpi:nbody-with-center.cpp user-defined.hpp  $(OBJS)
	$(PG_BUILD)
	$(MPICC) $(INC) $(MPICFLAGS) -DQUAD -o $@ nbody-with-center.cpp $(CLIBS)
nbody-with-center-quad:nbody-with-center.cpp user-defined.hpp  $(OBJS)
	$(PG_BUILD)
	$(CC) $(INC) $(CFLAGS) -pg -DQUAD -o $@ nbody-with-center.cpp $(CLIBS)

clean:
	rm -f *.o *~ nbody-with-center

distclean: clean
	$(PG_CLEAN)
	rm -f nbody-with-center
	rm -rf result


test: 
	# This command is only for FDPS developers.
	./test.py

export:
	rsync -avu $(SRCS) $(EXPORTDIR)
export-git:
	make export
	cd $(EXPORTDIR); git commit -a ; git push
# fdps-autotest-run (DO NOT CHANGE THIS LINE)
