#include "messages.hpp"

#if MPI_BUILD
    #include <mpi.h>
#endif
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.hpp"

/* An error function */
void error(Simulation& simulation, char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    fprintf(stderr, "+- ERROR (on proc %d): %s\n", simulation.me_, ws);
#if MPI_BUILD
    MPI_Finalize();
#endif
    exit(1);
}

/* And a warning function */
void warning(Simulation& simulation, char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    if (simulation.onRootProcessor()) {
        fprintf(stderr, "+- WARNING: %s\n", ws);
    }
}

/* Print some debug information to the screen */
void debugline(Simulation& simulation, int proc, char *message, ...) {
    va_list arg;
    static char ws[LINE_LENGTH];
    va_start(arg, message);
    vsprintf(ws, message, arg);
    va_end(arg);
    if (simulation.me_ == proc) {
        fprintf(stdout, "+- %s\n", ws);
    }
}

// void dumpToFiles(BUFFER *sendbuf, BUFFER *recvbuf, int bufsize) {

    // int i, f, nsr, *curbufsize, *lcbs;
// #if MPI_BUILD
    // MPI_Status mpistatus;
// #endif

    // [> Build an array of the current buffer size for broadcast <]
    // curbufsize = (int *)malloc(NProcessors*sizeof(int));
    // lcbs = (int *)malloc(NProcessors*sizeof(int));
    // for (i=0; i<NProcessors; ++i) {
        // lcbs[i] = 0;
    // }
    // lcbs[Me] = bufsize;
// #if MPI_BUILD
    // MPI_Allreduce(lcbs,curbufsize,NProcessors,MPI_INT,MPI_SUM,Universe);
// #else
    // curbufsize[Me] = lcbs[Me];
// #endif

// #if MPI_BUILD
    // [> Send the data to the root processor <]
    // if (Me != RootProc) {
        // MPI_Send(sendbuf,curbufsize[Me]*sizeof(BUFFER),MPI_CHAR,RootProc,0,Universe);
    // }
    // [> Receive the date from the daughter processors and write to file <]
    // else {
// #endif
        // [> Loop the processors <]
        // for (i=0; i<NProcessors; ++i) {
            // [> On the main processor we have to copy the data only, no send and receive <]
            // if (i==0) {
                // recvbuf = sendbuf;
            // }
// #if MPI_BUILD
            // [> For all other processors we need to receive the data <]
            // else {
                // MPI_Recv(recvbuf,curbufsize[i]*sizeof(BUFFER),MPI_CHAR,i,0,Universe,&mpistatus);
            // }
// #endif
            // [> Write data to file (only the root processor can do this) <]
            // [> PLEASE NOTE: DATA IS SENT IN STRIPED FORM, THEY ARE NOT ORDERED! <]
            // for (nsr=0; nsr<curbufsize[i]; ++nsr) {
                // f = recvbuf[nsr].iz;
                // [> The file buffer can be a gzip pipe or an ASCII file stream <]
                // fprintf(FStreams[f],"%d %d %d ",recvbuf[nsr].iz,recvbuf[nsr].ix,recvbuf[nsr].iy);
                // fprintf(FStreams[f],"%6.3f %6.3f %6.3f ",recvbuf[nsr].pos.x,recvbuf[nsr].pos.y,recvbuf[nsr].pos.z);
                // fprintf(FStreams[f],"%8.4f %8.4f %8.4f ",recvbuf[nsr].f.x,recvbuf[nsr].f.y,recvbuf[nsr].f.z);
                // fprintf(FStreams[f],"%6.3f %6.3f %6.3f ",recvbuf[nsr].d.x,recvbuf[nsr].d.y,recvbuf[nsr].d.z);
                // fprintf(FStreams[f],"%6.3f %8.4f ",recvbuf[nsr].dd,recvbuf[nsr].angle);
                // fprintf(FStreams[f],"%8.4f %d\n",recvbuf[nsr].e,recvbuf[nsr].n);
            // }
        // }
// #if MPI_BUILD
    // }
// #endif

    // [> Get rid of the buffer size broadcast arrays <]
    // free(curbufsize);
    // free(lcbs);

    // [> Go home <]
    // return;
// }
