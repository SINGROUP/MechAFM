#include "messages.hpp"

#if MPI_BUILD
    #include <mpi.h>
#endif
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.hpp"

/* An error function */
void error(char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    int process = 0;
#if MPI_BUILD
    MPI_Comm_rank(MPI_COMM_WORLD, &process);
#endif
    fprintf(stderr, "+- ERROR (on process %d): %s\n", process, ws);
#if MPI_BUILD
    MPI_Finalize();
#endif
    exit(1);
}

/* And a warning function */
void warning(char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    int process = 0;
#if MPI_BUILD
    MPI_Comm_rank(MPI_COMM_WORLD, &process);
#endif
    if (process == 0) {
        fprintf(stderr, "+- WARNING: %s\n", ws);
    }
}

/* Print a message to the user if we're the root.*/
void pretty_print(char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    int process = 0;
#if MPI_BUILD
    MPI_Comm_rank(MPI_COMM_WORLD, &process);
#endif
    if (process == 0) {
        fprintf(stdout, "+- %s\n", ws);
    }
}
