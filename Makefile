#------------------------------------------------------------------------------
# ROOT definitions
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lMathMore -lMinuit2

#------------------------------------------------------------------------------
# g++ options

CXX = g++
CXXFLAGS = -g -O -fPIC -Wall -Wno-deprecated
LD = g++

SOFLAGS = #-shared #-dynamiclib #-single_module #-undefined dynamic_lookup

CXXFLAGS     += $(ROOTCFLAGS)

LIBS         = ${ROOTGLIBS} ${ROOTLIBS}

#------------------------------------------------------------------------------
# Compilation

all : ./test/test

clean:
	rm *.o
	rm test*

objects= doRegress.o

./test/test: $(objects)
	   $(CXX) -o ./test/test $(objects) ${LIBS}  $(CXXFLAGS) ${SOFLAGS}

doRegress.o: doRegress.cpp base.h base_plot.h base_function.h simpleBP.h
