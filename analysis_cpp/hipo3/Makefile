ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

Clas12Tool := /home/sangbaek/C-BAND
cdir = $(shell pwd)

HIPOCFLAGS  := -I$(cdir) -I$(Clas12Tool)/Hipo -I$(Clas12Tool)/Banks -I$(Clas12Tool)/Physics
HIPOLIBS    := -L$(Clas12Tool)/lib -lhipo -lclas12banks -lclasphysics

LZ4LIBS     := -L/group/clas12/packages/lz4/lib -llz4
LZ4INCLUDES := -I/group/clas12/packages/lz4/lib

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

all: ep

ep: ep.o
	$(CXX) -o ep $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o ep *~

%.o: %.cpp
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
