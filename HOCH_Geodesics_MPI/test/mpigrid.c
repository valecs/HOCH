#include "mpigrid.h"

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>

#define BUF 128

/* A wrapper to elide the return value */
void MPI_Finalize_nr(void);

/*
  Parses command line arguments and returns a task ID based on the MPI
  rank. Expects an argument of the form -t n-m to sepcify tasks to be
  completed.  Inits MPI functions and arranges for MPI_Finalize to be
  called at exit.

  Idea: take a [-w] argument to redirect stdout and stderr to the
  specified working directories.
*/
int getTaskID(int argc, char ** argv){
  int n = 0, m = 0;

  int c = -1, s = 1;
  opterr = 0;  //since we do our own error-handling
  const char optstring[]="-t:";
  while (s && (-1 != (c = getopt(argc, argv, optstring)))){
    switch (c){
    case 't':
      fprintf(stderr, "MPIGRID: saw %c\n", c);
      /*
	Short-circuit evaluation ensures n & m are both set prior to
	comparison.
      */
      if (2 != sscanf(optarg, "%d-%d", &n, &m) ||
	  n > m){
	fprintf(stderr, "Bad range specification. Exiting!\n");
	exit(c);
      }
      
      /* Hide these options from future invocations of getopt(3). */
      /*
	I think what needs to happen here is for the arg(s) to be
	permuted to then *end* of argv.
      */
      argv[optind-1][0] = '\0';
      if (&optarg[0] != &argv[optind-1][2]){// if there was a space between args
	argv[optind-2][0] = '\0';
      }

      s = 0;
      break;
    default:
      /* don't want to interfere with others' options */
      fprintf(stderr, "MPIGRID: saw %c\n", c);
      break;
    }
  }

  // Reset optind for subsequent calls to getopt.
  optind = 1;
  
  MPI_Init(NULL, NULL);
  
  if (atexit(&MPI_Finalize_nr)){
    fprintf(stderr, "Failed to register exit handler! This is probably a problem.\n");
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  int taskID = rank + n;
  if (taskID > m){ // no need to run
    exit(0);
  }

  if (0 == rank){
    int nprocs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    if (nprocs < m - n + 1){
      fprintf(stderr, "Insufficient ranks for jobs requested;"
	      " the last job that will run is: %d\n", n + nprocs - 1);
    }
  }
  
  return taskID;
}

void MPI_Finalize_nr(void){
  MPI_Finalize();
}
