#include "mpigrid.h"

#include <mpi.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void MPI_Finalize_nr(void);

/*
  Parses command line arguments and returns a task ID based on the MPI
  rank. Expects an argument of the form -f n-m to sepcify tasks to be
  completed.  Inits MPI functions and arranges for MPI_Finalize to be
  called at exit.

  If a '-f' option is present, getTaskID will shift the contents of
  argv and change the count in argc to hide the elements from future
  invocations of getopt(3).

  If a '-f' option is present and the range cannot be parsed, the
  program will emit an error and terminate.

  If there are more tasks than ranks, rank 0 will emit a warning. If
  there are more ranks than tasks, the superflous ranks will be killed
  silently.
*/
int getTaskID(const char f, int * argc, char ** argv){
  int n = 0, m = 0; /* Assign tasks on [n,m] */

  /* prevent permuting argv via intiail '-' */
  char optstring[] = "-_:";
  /* Insert flag into option sting. */
  optstring[1] = f;
  
  int c = -1, optind_prev = optind;
  opterr = 0;  /* since we do our own error-handling */
  while (-1 != (c = getopt(*argc, argv, optstring))){
    if (f == c){
      /*
	Short-circuit evaluation ensures n & m are both set prior to
	comparison.
      */
      if (2 != sscanf(optarg, "%d-%d", &n, &m) ||
	  n > m){
	fprintf(stderr, "Bad range specification. Exiting!\n");
	exit(c);
      }

      /* Hide these elements from future invocations of getopt(3). */
      memmove(&argv[optind_prev],
	      &argv[optind],
	      (size_t) (*argc - optind) * sizeof(char *));
      *argc -= optind - optind_prev;

      /* Stop after this round. */
      break;
    }
    optind_prev = optind;
  }
  
  /* Reset optind for subsequent calls to getopt. */
  optind = 1;

  /* Setup the MPI enviroment */
  MPI_Init(NULL, NULL);
  
  if (atexit(&MPI_Finalize_nr)){
    fprintf(stderr, "Failed to register MPI exit handler! This is probably a problem.\n");
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

/*
  A wrapper around MPI_Finalize with the correct type signature for
  registration via atexit(3).
*/
void MPI_Finalize_nr(void){
  MPI_Finalize();
}
