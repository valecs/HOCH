#include "../src/kinds.h"
#include "../src/main.h"

#include <ctype.h>

#include <stdio.h>
#include <stdlib.h>
#include "mpigrid.h"
#include <getopt.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#define BUF 128

void msg(const char * s);
void getTaskFromFile(const char * fname);
void parseTask(char * out);
int testFileLoad(int argc, char ** argv);
void printArgs(int argc, char ** argv);

int manifestID = -1;

struct params p = {
  .manifestID = -1,
  .kind = nothing,
  .E0 = -1,
  .stratt_root = NULL,
};


int main(int argc, char ** argv){
  testFileLoad(argc, argv);
  
  return 0;
}

void printArgs(int argc, char ** argv){
  for (int i = 0; i < argc; i++){
    printf("%d: %s\n", i, argv[i]);
  }
}

void msg(const char * s){
  fprintf(stderr, "%s\n", s);
}


int testFileLoad(int argc, char ** argv){
  if (argc != 3){
    msg("Require 2 arguments.");
    return 1;
  }

  p.manifestID = atoi(argv[2]);
  fprintf(stderr, "Got: %d\n", p.manifestID);

  getTaskFromFile(argv[1]);

  printf("manifestID = %05d\n"
	 "kind = %s\n"
	 "energy = %05d\n",
	 p.manifestID, k2s(p.kind), p.E0);

  return 0;
}


/*
  Set a new manifestID, kind, and energy by reading the manifestID-th
  line from the file at path fname.  Clearly, this requires that
  manifestID already be set.
*/
void getTaskFromFile(const char * fname){
  
  FILE * f;
  f = fopen(fname, "r");
  if (!f){
    msg(fname);
    msg("FAILURE_READ_OPEN");
    exit(1);
  }

  char in[BUF];
  char out[BUF];
  int j = 0;
  int lines = p.manifestID;
  while (lines > 0 && 0 != fread(in, sizeof(char), BUF, f)){
    for(char * i = in; i < in + BUF ; i++){
      /*
	if *this* is our line,
	copy to output buffer
      */
      if (1 == lines){
	out[j++] = (*i);
      }

      if (j >= BUF){
	msg(fname);
	msg("Buffer error while reading task.");
	exit(1);
	break;
      }
      
      /*
	decrement lines on each '\n'
	if we hit 0, read out that line
      */
      if ('\n' == (*i) && !--lines){
	out[j-1] = '\0';
	parseTask(out);
	break;
      }
    }
  }

  fclose(f);
}


/*
  Helper-function for getTaskFromFile
*/
void parseTask(char * out){
  const char delim[] = " \t";
  char * tok;

  if (!(tok = strtok(out, delim))){return;}
  errno = 0;
  p.manifestID = (int) strtol(tok, (char **) NULL, 10);
  if (errno){
    msg("Cannot parse ID; exiting!\n");
    exit(1);
  }
	  
  if (!(tok = strtok(NULL, delim))){return;}
  p.kind = c2k(tok[0]);

  if (!(tok = strtok(NULL, delim))){return;}
  errno = 0;
  p.E0 = (int) strtol(tok, (char **) NULL, 10);
  if (errno){
    msg("Cannot parse energy; exiting!\n");
    exit(1);
  }
}

/*
  The memory returned is static and should not be freed
*/
const char * k2s(geodesicKind k){
  static char s[BUF];
  
  switch (k){
  case direct:
    sprintf(s, "direct");
    break;
  case roaming:
    sprintf(s, "roaming");
    break;
  case radical:
    sprintf(s, "radical");
    break;
  case nothing:
    sprintf(s, "nothing");
    break;
  default:
    sprintf(s, "error");
    break;
  }

  return s;
}

char k2c(geodesicKind k){
  switch (k){
  case direct:
    return 'D';
    break;
  case roaming:
    return 'R';
    break;
  case radical:
    return 'A';
    break;
  case nothing:
    return 'N';
    break;
  default:
    return '0';
    break;
  }
}


geodesicKind c2k(char c){
  switch (c){
  case 'D': /* exclicit fallthrough */
  case 'd':
    return direct;
    break;
  case 'R': /* exclicit fallthrough */
  case 'r':
    return roaming;
    break;
  case 'A': /* exclicit fallthrough */
  case 'a':
    return radical;
    break;
  case 'N': /* exclicit fallthrough */
  case 'n':
  default:
    return nothing;
    break;
  }
}
