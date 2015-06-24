#pragma once

#include "globals.hpp"
#include "simulation.hpp"

void error(char *message, ...);
void warning(char *message, ...);
void pretty_print(char *message, ...);
//void dumpToFiles(Simulation& simulation, BUFFER *sendbuf, BUFFER *recvbuf, int bufsize);
