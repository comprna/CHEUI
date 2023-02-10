EXEC=CHEUI
CC=g++
CFLAGS=-Wall -Wextra -O3 -fno-strict-aliasing -DLINUX_ -DOC_NEW_STYLE_INCLUDES -Wno-deprecated -I./PicklingTools/opencontainers_1_8_5/include -I./PicklingTools -pthread -D_REENTRANT

all: $(EXEC)

$(EXEC): main.cpp 
	$(CC) -o $(EXEC) $(CFLAGS) main.cpp kmer.hpp PicklingTools/libptools.so 

clean: 
	rm -f *.o
	rm -f $(EXEC)
	
.PHONY: all clean