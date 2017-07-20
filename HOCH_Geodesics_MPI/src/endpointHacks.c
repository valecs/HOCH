
/*
  some bit-rotted, nasty code for loading boundary conditions,
  hidden away in a separate translation unit.

  exports one function:
  int brittleLoad(gsl_vector * educt, gsl_vector * roamer, gsl_vector * product, int taskID, bool sampleRegion, bool sampleRadicalDissociation)
*/

#include "endpointHacks.h"
#include "kinds.h"
#include <linkedList.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h> // for INT_MAX

// from HOCHgeodesics.c
extern void path_free(llist * path);

// from main.c
extern void msg(const char * s);

static int  ext_gsl_vector_fscanfn(FILE * stream, gsl_vector * v, int n);
static int ext_gsl_vector_fscanf_safe(FILE * stream, gsl_vector *v);
static llist * fscanf_path(FILE * f, int skip, size_t dim);
static int fskipn(FILE * f, int n);

#define BUF 1024

/*
  manifest stored at:
    $STRATT_ROOT/data/HOCH/manifest/varel/E0.manifest
  with line structure:
    EDUCT PRODUCT_MOL PRODUCT_RAD
  roamers are colocated at:
    $STRATT_ROOT/data/HOCH/manifest/varel/E0.roamers
*/
void brittleLoad2(struct pointsHolder * points, struct params p){
  char fname[BUF];
  if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/manifest/varel/%05d.manifest",
	       p.stratt_root, p.E0) >= BUF ){
    msg("FAILURE_WRITE_TRUNCATE");
    exit(3);
  }

  FILE * f;
  f = fopen(fname, "r");
  if (!f){
    msg(fname);
    msg("FAILURE_READ_OPEN_MANIFEST");
    exit(3);
  }

  int ts_educt = 0, ts_prodMol = 0, ts_prodRad = 0;
  
  fskipn(f, (p.manifestID-1));
  if (3 != fscanf(f,"%05d %05d %05d\n", &ts_educt, &ts_prodMol, &ts_prodRad)){
    msg("FAILURE_PARSE_MANIFEST");
  }

  ////////////////// educt //////////////////
  if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/geometries/%05d-%05d.gamma",
	       p.stratt_root, p.E0, ts_educt) >= BUF ){
    msg("FAILURE_WRITE_TRUNCATE");
    exit(3);
  }
  
  f = fopen(fname,"r");
  if (!f){
    msg(fname);
    msg("FAILURE_READ_OPEN_EDUCT");
    exit(3);
  }
  // rely on position coming first in phase space representation.
  if (GSL_SUCCESS !=  gsl_vector_fscanf(f, points->A)){
    msg(fname);
    msg("FAILURE_PARSE_EDUCT");
    exit(3);
  }

  fclose(f);
  
  ////////////////// product //////////////////
  int ts_product = ts_prodMol;
  if (radical == p.kind){
    ts_product = ts_prodRad;
  }
  
  if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/trajectories/%05d/%05d.timeseries",
	       p.stratt_root, p.E0, ts_product) >= BUF ){
    msg("FAILURE_WRITE_TRUNCATE");
    exit(3);
  }

  f = fopen(fname,"r");
  if (!f){
    msg(fname);
    msg("FAILURE_READ_OPEN_PRODUCT");
    exit(3);
  }
  
  {
    llist * r = fscanf_path (f, 1, (size_t) 24);//since grabbing from phase space
    fclose(f);
    r = last(r);
    gsl_vector_view v = gsl_vector_subvector((gsl_vector *)(r->e), (size_t) 0, (size_t) 12);
    gsl_vector_memcpy(points->B, &v.vector);
    path_free(r);
  }

  ////////////////// roamer //////////////////
  if (roaming != p.kind){
    return;
  }

  // move B to C; roamer inserted at B
  points->C = points->B; points->B = NULL;

  if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/manifest/varel/%05d.roamers",
	       p.stratt_root, p.E0) >= BUF ){
    msg("FAILURE_WRITE_TRUNCATE");
    exit(3);
  }

  f = fopen(fname,"r");
  if (!f){
    msg(fname);
    msg("FAILURE_READ_OPEN_ROAMER");
    exit(3);
  }

  if (GSL_SUCCESS !=  ext_gsl_vector_fscanfn(f, points->B, p.manifestID)){
    msg(fname);
    msg("FAILURE_PARSE_ROAMER");
    exit(3);
  }

  fclose(f);    

}


// reads the nth one-indexed vector from stream. n=-1 specifies the last vector.
// can only input steams with up  to INT_MAX vectors, which on a x86_64 = 2^31 -1
// returns EXIT_FAILURE on fail, EXIT_SUCCESS otherwise
static int  ext_gsl_vector_fscanfn(FILE * stream, gsl_vector * v, int n){
  // -1 and the natural numbers are valid
  if ((n==0) || (n < -1))
    return EXIT_FAILURE;

  int limit;
  if (-1 == n)
    limit = INT_MAX;
  else
    limit = n;

  int retVal = EXIT_SUCCESS;

  //restore  gsl_set_error_handler(err_gsl);
  gsl_error_handler_t * err_gsl = gsl_set_error_handler_off();

  for (int i = 0; i<limit; i++){
    if (GSL_SUCCESS != ext_gsl_vector_fscanf_safe(stream,v)){
      if (-1 == n){
	break;
      }
      else{//happened too early!
	retVal = EXIT_FAILURE;
	break;
      }
    }
  }

  gsl_set_error_handler(err_gsl);

  return retVal;
}

// Identical to gsl_vector_fscanf except the vector v is only modified on success
static int ext_gsl_vector_fscanf_safe(FILE * stream, gsl_vector *v){
  gsl_vector * temp = gsl_vector_alloc(v->size);
  int retVal = gsl_vector_fscanf(stream, temp);

  if(GSL_SUCCESS == retVal)
    gsl_vector_memcpy(v,temp);

  gsl_vector_free(temp);
  return retVal;
}

///CRAP ALERT!!///

static int fskipn(FILE * f, int n){
  int count;
  int c;
  if (n <= 0){return 0;} //null case

  count = 0;
  while(count < n){
    c = fgetc(f);
    switch (c){
    case EOF:
      return -1;
      break;
    case (int)('\n'):
      count++;
      break;
    default:
      break;
    }
  }

  return count;
}

static llist * fscanf_path(FILE * f, int skip, size_t dim){
  gsl_vector * v, * t;
  int STATUS_V, STATUS_S;

  llist * l = ll_new(NULL);
  t = gsl_vector_alloc(dim);

  gsl_error_handler_t * err_gsl = gsl_set_error_handler_off(); //restore with gsl_set_error_handler(err_gsl);
  while(1){

    STATUS_S = fskipn(f, skip);

    STATUS_V = gsl_vector_fscanf(f,t);

    if ((STATUS_V != GSL_EFAILED) && (STATUS_S == ((int)skip) )){
      v = gsl_vector_alloc(dim);
      gsl_vector_memcpy(v,t);
      append(l,v);
      fskipn(f,1);
    }else{break;}

  }
  gsl_set_error_handler(err_gsl);
  
  gsl_vector_free(t);

  l = tail(l); //clips NULL element at begining;
  return l;
}
