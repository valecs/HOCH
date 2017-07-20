#include "HOCHgeodesics.h"
/* llist * getOptimizedHOCHGeodesic(
     const gsl_vector * A,
     const gsl_vector * B,
     double landscapeEnergy);
 */

#include "endpointHacks.h"
/*
   void brittleLoad2(
     struct pointsHolder * points,
     struct params p);
*/

#include "mpigrid.h"
#include "kinds.h"

#include <formaldehydeProperties.h>
#include <linkedList.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <getopt.h>

#define STRATT_ROOT "STRATT_ROOT"
#define HACM 219474.6313705
#define BUF 1024

struct params p = {
  .manifestID = -1,
  .kind = nothing,
  .E0 = -1,
  .stratt_root = NULL,
};

void usage(void);

struct pointsHolder loadEndPoints(void);
void setCM(struct pointsHolder points);
void writeRoute(llist * route);
void writePath(FILE * s, llist * p);
void * printVectorM(void * e);
void getTaskFromFile(const char * fname);

const char * k2s(geodesicKind k);
char k2c(geodesicKind k);
void msg(const char * s);
void * printVectorM(void * e);
/* only touched by printVectorM*/
static FILE * printVectorM_stream = NULL;

/*
  FIXME: Need a flag to also get ID from a file or a list.
*/
int main(int argc, char ** argv){
  p.stratt_root = getenv(STRATT_ROOT);
  p.manifestID = getTaskID(argc, argv);

  if (!p.stratt_root){
    fprintf(stderr, "Variable, " STRATT_ROOT ", not found; exiting!\n");
    exit(1);
  }

  {
    int c = 1;
    opterr = 0;  //since we do our own error-handling
    const char optstring[]="+e:hDRAF:";
    while (-1 != (c = getopt(argc, argv, optstring))){
      switch (c){
      case 'e':
	p.E0 = atoi(optarg);
	break;
      case 'D':
	p.kind = direct;
	break;
      case 'R':
	p.kind = roaming;
	break;
      case 'A':
	p.kind = radical;
	break;
      case 'F':
	getTaskFromFile(optarg);
	break;
      case 'h':
	usage();
	exit(0);
	break;
      case '?':
	msg("Bad argument; exiting!\n");
	exit(1);
	break;
      default:
	msg("Internal error; exiting!");
	exit(1);
	break;
      }
    }
  }
  
  // Validate parameters
  if (nothing == p.kind ||
      -1      == p.E0   ||
      0       >  p.manifestID){
    msg("Invalid parameters; exiting!");
    usage();
    exit(1);
  }

  const int E0 = p.E0;
  const geodesicKind kind = p.kind;
  
  struct pointsHolder points = loadEndPoints();

  setCM(points);
  
  llist * route = NULL;
  llist * AB = getOptimizedHOCHGeodesic(points.A , points.B , E0/HACM);

  if (!AB){
    msg("SUCCESS NR_AB");
    exit(0);
  }
  
  if (kind == direct ||
      kind == radical){
    route = AB;
    AB = NULL;
  }
  else{
    llist * BC = getOptimizedHOCHGeodesic(points.B , points.C , E0/HACM);
    if (!BC){
      msg("SUCCESS NR_BC");
      exit(2);
    }
    else{
      route = join(AB, tail(BC));
      AB = NULL; BC = NULL;
    }
  }

  writeRoute(route);  
  msg("SUCCESS");

  return 0;
}

void usage(void){
  fprintf(stderr, "HOCH_Geodesics_MPI -t n-m -e EL {-D | -R | -A} [-F TASKFILE]\n");
  fprintf(stderr, "  Computes a geodesic at landcape energy EL on the HOCH PES.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Parameters:\n");
  fprintf(stderr, "  EL : energy in wavenumbers.\n");
  fprintf(stderr, "  n-m: range of initial conditions to generate geodesics over; [n,m].\n");
  fprintf(stderr, "The final Flag Controls the kind of geodesic sought.\n");
  fprintf(stderr, "  -D Direct\n");
  fprintf(stderr, "  -R Roaming\n");
  fprintf(stderr, "  -A rAdical\n");
  fprintf(stderr, "The optional flag, -F causes jobs to be read from TASKFILE. -t is\n");
  fprintf(stderr, "  still required.\n");
  fprintf(stderr, "Requres the variable " STRATT_ROOT " be set.\n");
}

void msg(const char * s){
  fprintf(stderr, "%05d %c : %s\n", p.manifestID, k2c(p.kind), s);
}

void setCM(struct pointsHolder points){
  if (points.A){resetCenterOfMass(points.A);}
  if (points.B){resetCenterOfMass(points.B);}
  if (points.C){resetCenterOfMass(points.C);}
}

struct pointsHolder loadEndPoints(void){
  struct pointsHolder points = {
    gsl_vector_alloc(formaldehydeNumDimensions),
    gsl_vector_alloc(formaldehydeNumDimensions),
    gsl_vector_alloc(formaldehydeNumDimensions)
  };
  
  brittleLoad2(&points, p);
  return points;
}

/*
  Set a new manifestID by reading the manifestID-th line from the file at path fname.
  Clearly, this requires that manifestID already be set. 
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
	  p.manifestID = atoi(out);
	  break;
      }
    }
  }

  fclose(f);
}

void writeRoute(llist * route){
  char fname[BUF];
  if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/geodesics/varel/%05d/%05d.%s",
	       p.stratt_root, p.E0, p.manifestID, k2s(p.kind)) >= BUF ){
    msg("FAILURE_WRITE_TRUNCATE");
    exit(3);
  }

  FILE * f;
  f = fopen(fname, "w");
  if (!f){
    msg(fname);
    msg("FAILURE_WRITE_OPEN");
    exit(3);
  }

  printVectorM_stream = f;
  map(route, printVectorM);
  fclose(f);
}

/*
  prints a vector on printVectorM_stream; use like:
    printVectorM_stream = f;
    map(route, printVectorM);
*/
void * printVectorM(void * e){
  gsl_vector * v = (gsl_vector *)e;
  gsl_vector_fprintf(printVectorM_stream, v, "%#22.16g");
  return NULL;
}

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
