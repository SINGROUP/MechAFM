/* Load system headers */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <glob.h>

/* Some macro definitions */
#define LINE_LENGTH 512
#define NAME_LENGTH 128
#define ATOM_LENGTH 8
#define TOLERANCE 1e-10
#define TRUE 1
#define FALSE 0
#define NEGVAL -99999
#define PI 3.14159265358979323846
#define SIXTHRT2 1.12246204830937298142

/* Define vector types */
typedef struct vector {
  double x, y, z;
} VECTOR;
VECTOR NULL_vector = {0.0,0.0,0.0};
typedef struct ivector {
  int x, y, z;
} IVECTOR;
IVECTOR NULL_ivector = {0,0,0};


/**********************
 ** GLOBAL VARIABLES **
 **********************/

double Amplitude;
double KCantilever;
double Frequency0;
char FileBase[NAME_LENGTH];
int NFiles;
char **FileList;
double *Zpos;
VECTOR Delta;

/**********************
 ** HEADER FUNCTIONS **
 **********************/

/* The help function */
void help() {
  fprintf(stderr, "+-\n");
  fprintf(stderr, "+- Usage: convert-forces-to-deltaf <options> filebase <options>\n");
  fprintf(stderr, "+-        -a ##        amplitude (in Angstrom)\n");
  fprintf(stderr, "+-        -k ##        cantilever stiffness (in kcal/mol/Angstrom/Angstrom)\n");
  fprintf(stderr, "+-        -f ##        cantilever oscillating frequency (in Hz)\n");
  fprintf(stderr, "+-\n");
  return;
}

/* An error function */
void error(char *message, ...) { 
  va_list arg;
  static char ws[LINE_LENGTH];
  va_start(arg, message);
  vsprintf(ws, message, arg);
  va_end(arg);
  fprintf(stderr, "+- ERROR: %s\n", ws);
  help();
  exit(1);
}

/* And a warning function */
void warning(char *message, ...) { 
  va_list arg;
  static char ws[LINE_LENGTH];
  va_start(arg, message);
  vsprintf(ws, message, arg);
  va_end(arg);
  fprintf(stderr, "+- WARNING: %s\n", ws);
}

/* Print some debug information to the screen */
void debugline(char *message, ...) {
  va_list arg;
  static char ws[LINE_LENGTH];
  va_start(arg, message);
  vsprintf(ws, message, arg);
  va_end(arg);
  fprintf(stdout, "+- %s\n", ws);
}

/* Convert a string to uppercase */
char *strupp(char *string) {
  char *convert;
  convert = string;
  do { *convert = toupper((unsigned char)*convert); } while (*convert++);
  return string;
}

/* Convert a string to lowercase */
char *strlow(char *string) {
  char *convert;
  convert = string;
  do { *convert = tolower((unsigned char)*convert); } while (*convert++);
  return string;
}
  
/* Filter out comment and empty lines in a file */
int checkForComments(char *line) {
  int moveon = FALSE;
  if (line[0] == '#') { moveon = TRUE; }
  else if (line[0] == '%') { moveon = TRUE; }
  else if (line[0] == '\n') { moveon = TRUE; }
  return moveon;
}

/* Check if a value is an integer */
int isint(char *str) {
  int integer = TRUE;
  int n = strlen(str);
  int i;
  for (i=0; i<n; ++i) {
    if (isdigit(str[i]) || str[i] == '-' || str[i] == '+') continue;
    integer = FALSE;
  }
  return integer;
}

/* Get the distance from the file name */
double getNumber(const char *str) {
  while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1))))) {
    str++;
  }
  return strtod(str,NULL);
}

/* Get file extension */
const char *getFileExt(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

/**************************
 ** FILE INPUT FUNCTIONS **
 **************************/

/* Read stuff from the command line */
void parseCommandLine(int argc, char *argv[]) {
  double amplitude = -1.0;
  double kcantilever = -1.0;
  double frequency0 = -1.0;
  char filename[NAME_LENGTH];
  int index, c, k;

  opterr = 0;

  /* Check the options */
  while ((c = getopt (argc, argv, "a:k:f:")) != -1) {
    switch (c) {
      case 'a':
        amplitude = atof(optarg);
        break;
      case 'k':
	kcantilever = atof(optarg);
        break;
      case 'f':
        frequency0 = atof(optarg);
        break;
      case '?':
        if ((optopt == 'a') || (optopt == 'k') || (optopt == 'f'))
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        return;
      default:
        abort ();
    }
  }
  
  /* Parse the remaining arguments */
  for (index = optind, k = 0; index < argc; index++) {
    if (k == 0) { sprintf(filename,"%s",argv[index]); }
    else { error("Invalid argument %s\n",argv[index]); }
    k++;
  }

  /* Store them to the global variables */
  Amplitude = amplitude;
  KCantilever = kcantilever;
  Frequency0 = frequency0;
  sprintf(FileBase,"%s",filename);

  return;
}

void printOptions() {
  fprintf(stderr, "+- Conversion of forces to frequency shifts (part of the MechAFM)\n");
  fprintf(stderr, "+-\n");
  fprintf(stderr, "+- Amplitude (Angstrom)                              = %12.3f\n",Amplitude);
  fprintf(stderr, "+- Cantilever stiffness (kcal/mol/Angstrom/Angstrom) = %12.3f\n",KCantilever);
  fprintf(stderr, "+- Cantilever oscillating frequency (Hz)             = %12.3f\n",Frequency0);
  fprintf(stderr, "+- Base name of force files to be processed          = %12s\n",FileBase);
  fprintf(stderr, "+-\n");
  return;
}

/* A file pointer handling the gzip format */
FILE *openFile(const char *fname) {
  FILE *fp;
  char cmd[NAME_LENGTH];
  if (strcmp(getFileExt(fname),"gz")==0) {
    sprintf(cmd,"gunzip -c %s",fname);
    fp = popen(cmd,"r");
  }
  else {
    fp = fopen(fname,"r");
  }
  if (fp==NULL) { error("The file %s does not exist!",fname); }
  return fp;
}

/* Check our overall sanity */
void assessUniverse() {

  glob_t forcefiles;
  size_t cnt, length;
  char **p;
  char wildcard[NAME_LENGTH], cmd[NAME_LENGTH];
  char line[LINE_LENGTH],  dump[LINE_LENGTH];
  int i, nlines, di;
  FILE *fp;
  int ix, iy, oldix, oldiy;
  double x, y, oldx, oldy;


  /* Retrieve the names of all force files */
  sprintf(wildcard,"%s*",FileBase);
  glob(wildcard, GLOB_NOCHECK, 0, &forcefiles);
  NFiles = 0;
  for (p=forcefiles.gl_pathv, cnt=forcefiles.gl_pathc; cnt; p++, cnt--) { NFiles++; }
  FileList = (char **)malloc(NFiles*sizeof(char*));
  for (i=0; i<NFiles; ++i) { FileList[i] = (char *)malloc(NAME_LENGTH*sizeof(char)); }
  for (p=forcefiles.gl_pathv, cnt=forcefiles.gl_pathc, i=0; cnt; p++, cnt--, ++i) { 
    sprintf(FileList[i],"%s",*p);
  }

  /* Get an array of z values (from the file names) */
  Zpos = (double *)malloc(NFiles*sizeof(double));
  for (i=0; i<NFiles; ++i) { Zpos[i] = fabs(getNumber(FileList[i])); }
  Delta.z = Zpos[1]-Zpos[0];

  /* Get resolution information from the first file */
  fp = openFile(FileList[0]);
  while (fgets(line, LINE_LENGTH, fp)!=NULL) { nlines++; }
  fclose(fp);
  fp = openFile(FileList[0]);
  Delta.x = -1.0;
  Delta.y = -1.0;
  i = 0;
  while (fgets(line, LINE_LENGTH, fp)!=NULL) {
    sscanf(line,"%d %d %d %lf %lf %s",&di,&ix,&iy,&x,&y,dump);
    if (i==0) { oldx = x; oldy = y; oldix = ix; oldiy = iy; i++; continue; }
    if (ix == (oldix+1)) { Delta.x = x-oldx; }    
    if (iy == (oldiy+1)) { Delta.y = y-oldy; }
    i++;
    if ((Delta.x>0) && (Delta.y>0)) { break; }
  }

  return;
}

void readFiles() {
  return;
}

void writeData() {
  return;
}


/**********************
 ** THE MAIN ROUTINE **
 **********************/

int main(int argc, char *argv[]) {

  /* Read the command line */
  parseCommandLine(argc,argv);
  printOptions();
  
  /* Read the files */
  assessUniverse();
  readFiles();

  /* Store the data */
  writeData();

  /* Done */
  return 0;
}
