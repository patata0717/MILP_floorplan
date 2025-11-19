CC     = g++
CFLAGS = -std=c++11

OBJS   = main.o datastructure.o

# Default target: build both hw3 and visualize
all: hw3 visualize

# -------------------------------------------------------------------
# hw3 binary (your ILP generator)
# -------------------------------------------------------------------
hw3: $(OBJS)
	$(CC) $(OBJS) -o ../bin/hw3

main.o: main.cpp datastructure.hpp
	$(CC) -c main.cpp $(CFLAGS)

datastructure.o: datastructure.cpp datastructure.hpp
	$(CC) -c datastructure.cpp $(CFLAGS)

# -------------------------------------------------------------------
# visualize binary (SVG floorplan viewer)
#   Source:  ../visualize/visualize.cpp
#   Output:  ../bin/visualize
# -------------------------------------------------------------------
VIS_SRC = ../visualize/visualize.cpp
VIS_BIN = ../bin/visualize

visualize: $(VIS_BIN)

$(VIS_BIN): $(VIS_SRC)
	$(CC) $(CFLAGS) $< -o $@

# -------------------------------------------------------------------
# Cleanup
# -------------------------------------------------------------------
clean:
	rm -f *.o ../bin/hw3 ../bin/visualize
