CC = mpicc
NVCC = nvcc

SDKPATH := ~/NVIDIA_GPU_Computing_SDK
CUDAPATH := /usr/local/cuda

TARGET = smatcher

OBJS = kmp.o bm.o ac.o sh.o sbom.o wu.o sog8.o main.o helper.o cuda_ac.o cuda_sh.o cuda_sbom.o cuda_wm.o cuda_sog.o

#CPPFLAGS=-Wall -Wno-pointer-sign -O0 -g -funroll-loops -pg
CPPFLAGS=-Wall -Wno-pointer-sign -O2 -funroll-loops
NVCCFLAGS= -O2 -I$(CUDAPATH)/include -I$(SDKPATH)/shared/inc -I$(SDKPATH)/C/common/inc --ptxas-options=-v -arch=compute_12 -code=sm_12,compute_12

#LDFLAGS=-L$(CUDAPATH)/lib -L$(SDKPATH)/shared/lib/linux -L$(SDKPATH)/C/lib -lcuda -lcudart -lmpich
LDFLAGS=-L$(CUDAPATH)/lib -lcuda -lcudart -lmpich

all: $(TARGET)

$(TARGET): $(OBJS) $(SEQUENTIAL-OBJS)
	$(CC) $(CPPFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)
	
main.o: main.c
	$(CC) $(CPPFLAGS) -c main.c

kmp.o: kmp/kmp.c
	$(CC) $(CPPFLAGS) -c kmp/kmp.c
	
bm.o: bm/bm.c
	$(CC) $(CPPFLAGS) -c bm/bm.c
	
ac.o: ac/ac.c
	$(CC) $(CPPFLAGS) -c ac/ac.c
	
sh.o: sh/sh.c
	$(CC) $(CPPFLAGS) -c sh/sh.c
	
sbom.o: sbom/sbom.c
	$(CC) $(CPPFLAGS) -c sbom/sbom.c
	
wu.o: wu/wu.c
	$(CC) $(CPPFLAGS) -c wu/wu.c
	
sog8.o: sog/sog8.c
	$(CC) $(CPPFLAGS) -c sog/sog8.c

helper.o: ../helper.c
	$(CC) $(CPPFLAGS) -c ../helper.c

cuda_ac.o: cuda/cuda_ac.cu
	$(NVCC) $(NVCCFLAGS) -c cuda/cuda_ac.cu
	
cuda_sh.o: cuda/cuda_sh.cu
	$(NVCC) $(NVCCFLAGS) -c cuda/cuda_sh.cu
	
cuda_sbom.o: cuda/cuda_sbom.cu
	$(NVCC) $(NVCCFLAGS) -c cuda/cuda_sbom.cu
	
cuda_wm.o: cuda/cuda_wm.cu
	$(NVCC) $(NVCCFLAGS) -c cuda/cuda_wm.cu
	
cuda_sog.o: cuda/cuda_sog.cu
	$(NVCC) $(NVCCFLAGS) -c cuda/cuda_sog.cu

clean:
	rm -f *.o *.d $(TARGET) core
