CPP := g++

COMMON_FLAGS := -std=c++17 -DZI_USE_OPENMP -fopenmp -Wall -fPIC -O2 -march=native -mtune=native -funroll-loops -DNDEBUG 
LD_FLAGS := -Wl,-rpath,$(LOCAL_PATH)/armadillo/lib -L$(LOCAL_PATH)/armadillo/lib -lcblas -lblas -llapack -larmadillo
LD_FLAGS += -I/usr/local/include/ -I/usr/local/include/bamtools -I/usr/include/bamtools -I/home/pratth/.bamtools/usr/local/include -I/home/pratth/.bamtools/usr/local/include/bamtools -I/usr/include/bamtools -lbamtools -lz

OBJ_DIR := build
HEADERS := HEADERS = $(call rwildcard, src/, *.h) $(call rwildcard, src/, *.hpp)

all: $(OBJ_DIR) bin/mnasesignal
OBJ := build/bam.o build/gaussian.o

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p bin

build/bam.o: src/lib/bam.cpp
	$(CPP) $(LD_FLAGS) -fpermissive $(COMMON_OPT) -c -o $@ $<

build/gaussian.o: src/lib/gaussian.cpp
	$(CPP) $(LD_FLAGS) -fpermissive $(COMMON_OPT) -c -o $@ $<

bin/mnasesignal: $(OBJ) src/main.cpp
	$(CPP) -o $@ $^ $(COMMON_OPT) $(LD_FLAGS)
