#pragma once

#include "globals.hpp"

void error(char *message, ...);
void warning(char *message, ...);
void debugline(int proc, char *message, ...);
void dumpToFiles(BUFFER *sendbuf, BUFFER *recvbuf, int bufsize);
