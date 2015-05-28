## Parameters and files
MEXEC = mechafm-mpi
SEXEC = mechafm-serial
CFILES = mechafm-mpi.c

## Compiler
# On local machine
MCC = /usr/bin/mpicc
SCC = /usr/bin/gcc

## Local directory tree ##
INCDIR = .

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

############################################
## Actual make code below (do not change) ##
############################################

## Make the executable (MPI) ##
$(MEXEC): $(FILES)
    $(MCC) $(FULLFLAG) $(MPI_INC) $(MPI_PATH) $^ $(MATHFLAG) $(MPI_LIB) -o $(MEXEC)
    mkdir -p bin
    mv $(MEXEC) bin

## Make the executable (serial) ##
$(SEXEC): $(FILES)
    $(SCC) $(FULLFLAG) $(SERIAL) $^ $(MATHFLAG) -o $(SEXEC)
    mkdir -p bin
    mv $(SEXEC) bin

## Make clean ##
clean:
    rm -rf bin/$(MEXEC) bin/$(SEXEC) *~

## Make all ##
all: $(MEXEC) $(SEXEC)
