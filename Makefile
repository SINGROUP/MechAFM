## Parameters and files
MEXEC = mechafm-mpi
SEXEC = mechafm-serial
CFILES = src/mechafm-mpi.cpp src/messages.cpp src/messages.hpp

## Compiler
# On local machine
MCC = /usr/bin/mpic++
SCC = /usr/bin/c++

## Local directory tree ##
INCDIR = src/

SRCDIR = src/
BUILDDIR = build/
BINDIR = bin/
VPATH = $(SRCDIR) $(BUILDDIR) $(BINDIR)

## Flag settings ##
SERIAL   = -DSERIAL
OPTIM    = -O3 -fomit-frame-pointer
MATHFLAG = -lm
WARNFLAG = -Wshadow -Wno-format-zero-length -Wno-write-strings
FULLFLAG = $(OPTIM) $(WARNFLAG) -I$(INCDIR) -std=c++11

## Parallel thingies
MPI_INC  = -I/usr/lib/openmpi/include/
MPI_PATH = -L/usr/lib/openmpi/lib/
MPI_LIB  = #-lmpich

## Reshuffle all files ##
FILES = $(CFILES)
objects = main.o messages.o utility.o parse.o vectors.o system.o simulation.o

############################################
## Actual make code below (do not change) ##
############################################

## Make the executable (MPI) ##
$(MEXEC): $(objects)
	$(MCC) $(FULLFLAG) $^ $(MATHFLAG) $(MPI_LIB) -o $@
	mkdir -p $(BINDIR)
	mv $@ $(BINDIR)
	mkdir -p $(BUILDDIR)
	mv -u *.o $(BUILDDIR)

## Make the executable (serial) ##
$(SEXEC): $(objects)
	$(SCC) $(FULLFLAG) $(SERIAL) $^ $(MATHFLAG) -o $@
	mkdir -p $(BINDIR)
	mv $@ $(BINDIR)
	mkdir -p $(BUILDDIR)
	mv -u *.o $(BUILDDIR)

main.o: mechafm-mpi.cpp messages.hpp globals.hpp parse.hpp simulation.hpp system.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
simulation.o: simulation.cpp simulation.hpp globals.hpp system.hpp vectors.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
system.o: system.cpp system.hpp globals.hpp vectors.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
messages.o: messages.cpp messages.hpp globals.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
utility.o: utility.cpp utility.hpp globals.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
parse.o: parse.cpp parse.hpp globals.hpp messages.hpp utility.hpp vectors.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
vectors.o: vectors.cpp vectors.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@

## Make clean ##
clean:
	rm -f *.o
	rm -rf $(BINDIR)$(MEXEC) $(BINDIR)$(SEXEC)
	rm -rf $(BUILDDIR)

## Make all ##
all: $(MEXEC) $(SEXEC)
