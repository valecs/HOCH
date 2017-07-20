#include "pathHelpers.h"

#include <linkedList.h>

#include <gsl/gsl_blas.h>

void print_vect(gsl_vector * v){
  size_t dim = v->size; //cast from size_t
  for (size_t i = 0; i < (dim - 1); i ++){
    printf("%#16.14g\t",gsl_vector_get(v, i));
  }
  printf("%#16.14g\n",gsl_vector_get(v, dim-1));
}

void print_matrix(gsl_matrix * M){
  for   (size_t i = 0; i < M -> size1; i++){
    for (size_t j = 0; j < M -> size2; j++){
      printf("%#8.8g\t",gsl_matrix_get(M,i,j));
    }
    printf("\n");
  }
}

void print_path(llist * path){
  path = first(path);
  while (path->succ != NULL){
    print_vect(path->e);
    path = path->succ;
  }
  print_vect(path->e);
}

void fprintf_path (FILE * stream, llist * path, const char * format){
  path = first(path);
  if (path->e != NULL){
    do{
      gsl_vector_fprintf (stream, (gsl_vector *)(path->e), format);
    }while((path->succ != NULL) && (path = path->succ)); //order important: shortcircuit
  }
}

void print_lengths(llist * path){
  static size_t dim = 0; // assume we have no zero dimension spaces
  static gsl_vector * dummy   = NULL;

  if (dim   == 0)    {dim   = ((gsl_vector *) path->e)->size;} //determine dimensionality of space
  if (dummy == NULL) {dummy = gsl_vector_alloc (dim);}

  path = first(path);
  while (path->succ != NULL){
    gsl_vector_memcpy(dummy,path->e);
    gsl_vector_sub(dummy,path->succ->e);
    printf("%#16.14g\n",gsl_blas_dnrm2(dummy));
    path = path->succ;
  }
}

double path_length(llist * path){
  static size_t dim = 0; // assume we have no zero dimension spaces
  static gsl_vector * dummy   = NULL;
  double total_length = 0.0;

  if (dim   == 0)    {dim   = ((gsl_vector *) path->e)->size;} //determine dimensionality of space
  if (dummy == NULL) {dummy = gsl_vector_alloc (dim);}

  path = first(path);
  while (path->succ != NULL){
    gsl_vector_memcpy(dummy,path->e);
    gsl_vector_sub(dummy,path->succ->e);
    total_length += gsl_blas_dnrm2(dummy);
    path = path->succ;
  }
  return total_length;
}

//frees memory associated with a path
void free_path(llist * path){
  path = first(path);
  while (path->succ != NULL){
    gsl_vector_free(path->e); //free vector element
    path = path->succ;        //free container
    free(path->pred);
  }
  gsl_vector_free(path->e);
  free(path);
}
