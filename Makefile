## Parameters and files
MEXEC = mpi
SEXEC = serial

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
DEBUG    = -g
OPTIM    = -O3 -fomit-frame-pointer
MATHFLAG = -lm
WARNFLAG = -Wall -Wextra -Wshadow -Wno-format-zero-length -Wno-write-strings
FULLFLAG = $(DEBUG) $(OPTIM) $(WARNFLAG) -I$(INCDIR) -std=c++11

## Parallel thingies
MPI_BUILD = -D MPI_BUILD
MPI_INC  = -I/usr/lib/openmpi/include/
MPI_PATH = -L/usr/lib/openmpi/lib/
MPI_LIB  = #-lmpich
MPI_FLAGS = $(MPI_INC) $(MPI_PATH) $(MPI_LIB) $(MPI_BUILD)


mpi_objects = main-mpi.o messages-mpi.o simulation-mpi.o 
serial_objects = main-serial.o messages-serial.o simulation-serial.o 
shared_objects = parse.o system.o utility.o interactions.o minimiser.o integrators.o

############################################
## Actual make code below (do not change) ##
############################################

## Make the executable (MPI) ##
$(MEXEC): $(mpi_objects) $(shared_objects)
	$(MCC) $(FULLFLAG) $^ $(MATHFLAG) $(MPI_FLAGS) -o mechafm-$@
	mkdir -p $(BINDIR)
	mv mechafm-$@ $(BINDIR)
	mkdir -p $(BUILDDIR)
	mv -f *.o $(BUILDDIR)

## Make the executable (serial) ##
$(SEXEC): $(serial_objects) $(shared_objects)
	$(SCC) $(FULLFLAG) $^ $(MATHFLAG) -o mechafm-$@
	mkdir -p $(BINDIR)
	mv mechafm-$@ $(BINDIR)
	mkdir -p $(BUILDDIR)
	mv -f *.o $(BUILDDIR)

main-mpi.o: mechafm.cpp messages.hpp globals.hpp parse.hpp simulation.hpp system.hpp
	$(MCC) -c $(FULLFLAG) $< $(MATHFLAG) $(MPI_FLAGS) -o $@
main-serial.o: mechafm.cpp messages.hpp globals.hpp parse.hpp simulation.hpp system.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@
messages-mpi.o: messages.cpp messages.hpp globals.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_FLAGS) $< $(MATHFLAG) -o $@
messages-serial.o: messages.cpp messages.hpp globals.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@
simulation-mpi.o: simulation.cpp simulation.hpp globals.hpp system.hpp vectors.hpp minimiser.hpp integrators.hpp
	$(MCC) -c $(FULLFLAG) $(MPI_FLAGS) $< $(MATHFLAG) -o $@
simulation-serial.o: simulation.cpp simulation.hpp globals.hpp system.hpp vectors.hpp minimiser.hpp integrators.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@

system.o: system.cpp system.hpp globals.hpp vectors.hpp interactions.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@
utility.o: utility.cpp utility.hpp globals.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@
parse.o: parse.cpp parse.hpp globals.hpp messages.hpp utility.hpp vectors.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@
interactions.o: interactions.cpp interactions.hpp vectors.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@
minimiser.o: minimiser.cpp minimiser.hpp integrators.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@
integrators.o: integrators.cpp integrators.hpp
	$(SCC) -c $(FULLFLAG) $(MATHFLAG) $< -o $@

## Make clean ##
clean:
	rm -f *.o
	rm -f $(BINDIR)*
	rm -f $(BUILDDIR)*

## Make all ##
all: $(MEXEC) $(SEXEC)
