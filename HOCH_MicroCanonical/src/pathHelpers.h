#ifndef DISPLAYHELPERS_H
#define DISPLAYHELPERS_H

#include <linkedList.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

void print_vect(gsl_vector * v);
void print_matrix(gsl_matrix * M);
void print_path(llist * path);
void fprintf_path (FILE * stream, llist * path, const char * format);
void print_lengths(llist * path);
void free_path(llist * path);
double path_length(llist * path);

#endif
