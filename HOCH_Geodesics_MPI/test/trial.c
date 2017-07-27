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
int testFileLoad(int argc, char ** argv);
void printArgs(int argc, char ** argv);
void processOptions(int argc, char ** argv);

int manifestID = -1;

int main(int argc, char ** argv){
  printArgs(argc, argv);
  
  manifestID = getTaskID(argc, argv);
  printf("---Got manifestID = %d---\n", manifestID);
  printArgs(argc, argv);
  
  msg("---Begin option reporting---");
  processOptions(argc, argv);
  msg("---Final option run---");
  printArgs(argc, argv);
  return 0;
}

void processOptions(int argc, char ** argv){
  char loadFile[BUF] = {0};
  int c = 1;
  opterr = 0;  //since we do our own error-handling
  const char optstring[]="+e:hF:G";
  while (-1 != (c = getopt(argc, argv, optstring))){
    switch (c){
    case 'e':
      {
	errno = 0;
	int E = strtol(optarg, (char **) (NULL), 10);
	if (errno){
	  msg("Failed to read read energy");
	}
	else{
	  fprintf(stderr, "E=%d\n", E);
	}
      }
      break;
    case 'F':
      strncpy(loadFile, optarg, BUF);
      if ('\0' != loadFile[BUF-1]){
	msg(optarg);
	msg("File argument longer than BUF; things will break!");
      }
      else{
	msg(loadFile);
	msg("Read loadname.");
      }
      break;
    case 'G':
      {
	char * envvar = getenv("SGE_TASK_ID");
	if (!envvar){
	  msg("Variable, SGE_TASK_ID, not found; nothing to do.");
	  break;
	}
	  
	errno = 0;
	manifestID = strtol(envvar, (char **) (NULL), 10);
	if (errno){
	  msg("Cannot read SGE_TASK_ID.");
	}
	else{
	  fprintf(stderr, "Got manifestID=%d\n", manifestID);
	}
      }
      break;
    case 'h':
      msg("Help branch.");
      break;
    case '?':
      msg("Hit ? branch; bad argument.");
      break;
    default:
      fprintf(stderr, "Hit default branch; internal error. c = %o\n", c);
      break;
    }
  }
  if (loadFile[0]){
    /*
      getTaskFromFile requires that p.manifestID be set so wait
      until other options are processed.
    */
    getTaskFromFile(loadFile);
    fprintf(stderr, "Loaded from file %s and got %d.\n", loadFile, manifestID);
  }  
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

  manifestID = atoi(argv[2]);
  fprintf(stderr, "Got: %d\n", manifestID);

  getTaskFromFile(argv[1]);

  printf("manifestID = %05d\n", manifestID);

  return 0;
}

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
  int lines = manifestID;
  while (0 != fread(in, sizeof(char), BUF, f) && lines > 0){
    for(char * i = &in[0]; i < &in[BUF]; i++){
      if (1 == lines){
	out[j++] = (*i);
      }

      if (j >= BUF){
	msg("length error");
	return;
      }

      if ('\n' == (*i) && !--lines){
	  out[j-1] = '\0';
	  manifestID = atoi(out);
	  break;
      }
    }
  }

  fclose(f);
}
