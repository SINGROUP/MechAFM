#include "utility.hpp"

#include <ctype.h>
#include <string.h>

#include "globals.hpp"

char* strupp(char* string) {
    char *convert;
    convert = string;
    do {
        *convert = toupper((unsigned char) *convert);
    } while (*convert++);
    return string;
}

char* strlow(char* string) {
    char *convert;
    convert = string;
    do {
        *convert = tolower((unsigned char)*convert);
    } while (*convert++);
    return string;
}

int isint(char* str) {
    bool integer = true;
    int n = strlen(str);
    int i;
    for (i=0; i<n; ++i) {
        if (isdigit(str[i]) || str[i] == '-' || str[i] == '+'){
            continue;
        }
        integer = false;
    }
    return integer;
}
