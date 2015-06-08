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
FULLFLAG = $(OPTIM) $(WARNFLAG) -I$(INCDIR)

## Parallel thingies
MPI_INC  = -I/usr/lib/openmpi/include/
MPI_PATH = -L/usr/lib/openmpi/lib/
MPI_LIB  = #-lmpich

## Reshuffle all files ##
FILES = $(CFILES)
objects = main.o messages.o utility.o parse.o physics.o grid.o flexible.o vector.o

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
$(SEXEC): $(FILES)
	$(SCC) $(FULLFLAG) $(SERIAL) $^ $(MATHFLAG) -o $(SEXEC)
	mkdir -p $(BINDIR)
	mv $(SEXEC) $(BINDIR)

main.o: mechafm-mpi.cpp messages.hpp globals.hpp utility.hpp parse.hpp physics.hpp \
		grid.hpp flexible.hpp simulation.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
simulation.o: simulation.cpp simulation.hpp 
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
system.o: system.cpp system.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
messages.o: messages.cpp messages.hpp globals.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
utility.o: utility.cpp utility.hpp globals.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
parse.o: parse.cpp parse.hpp globals.hpp messages.hpp utility.hpp physics.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
physics.o: physics.cpp physics.hpp  globals.hpp grid.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
grid.o: grid.cpp grid.hpp globals.hpp messages.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
flexible.o: flexible.cpp flexible.hpp globals.hpp messages.hpp parse.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@
vector.o: vector.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $< $(MATHFLAG) $(MPI_LIB) -o $@

## Make clean ##
clean:
	rm -f *.o
	rm -rf $(BINDIR)$(MEXEC) $(BINDIR)$(SEXEC)
	rm -rf $(BUILDDIR)

## Make all ##
all: $(MEXEC) $(SEXEC)
