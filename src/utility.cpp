#include "utility.hpp"

#include <string.h>

#include "globals.hpp"

/* Convert a string to uppercase */
char *strupp(char *string) {
    char *convert;
    convert = string;
    do {
        *convert = toupper((unsigned char)*convert);
    } while (*convert++);
    return string;
}

/* Convert a string to lowercase */
char *strlow(char *string) {
    char *convert;
    convert = string;
    do {
        *convert = tolower((unsigned char)*convert);
    } while (*convert++);
    return string;
}

/* Check if a value is an integer */
int isint(char *str) {
    int integer = TRUE;
    int n = strlen(str);
    int i;
    for (i=0; i<n; ++i) {
        if (isdigit(str[i]) || str[i] == '-' || str[i] == '+'){
            continue;
        }
        integer = FALSE;
    }
    return integer;
}
