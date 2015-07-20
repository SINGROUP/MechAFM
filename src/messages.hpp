#pragma once

#include "globals.hpp"

// Prints an error message to the user and closes the program
void error(char* message, ...);
// Root process prints a warning message to the user
void warning(char* message, ...);
// Root process prints a message to the user
void pretty_print(char* message, ...);
