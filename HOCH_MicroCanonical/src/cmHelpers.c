#include "cmHelpers.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

#include <assert.h>

double det_3x3(gsl_matrix * M);

double kinetic_energy(gsl_vector * v, gsl_vector * m){
  assert(v->size==m->size);
  size_t dim = v->size;
  double e = 0;

  for (size_t i = 0; i < dim; i++){
    e += 0.5*(
	      gsl_vector_get(m,i) * 
	      gsl_pow_2(gsl_vector_get(v,i))
	      );
  }

  return e;
}
void get_inertia_tensor(gsl_vector * r, gsl_vector * m, gsl_matrix * I){
  const size_t space = 3; // gives the dimensionality of the space
  size_t dim = r->size;
  assert ((dim == r->size) && (dim == m->size));
  assert ((I->size1 == space) && (I->size2 == space));
  assert (dim % space == 0);
  size_t centers = dim / space;

  gsl_matrix_set_zero(I);

  gsl_vector * cm = gsl_vector_alloc(space);
  center_of_mass(r, m, cm);

  gsl_vector_const_view mass = gsl_vector_const_subvector_with_stride(m, (size_t) 0, space, centers);
  gsl_vector_view ri;

  //moves raw values from r to centers_mat
  gsl_matrix * centers_mat = gsl_matrix_calloc(centers, space);
  for (size_t row = 0; row < centers; row++){
    for (size_t col = 0; col < space; col++){
      gsl_matrix_set(centers_mat, row, col
		     ,gsl_vector_get(r, (row * space) + col));
    }
  }

  // shift to CM coors and construct mag^2
  gsl_vector * mag2 = gsl_vector_calloc(centers);
  for (size_t center = 0; center < centers; center++){
    ri = gsl_matrix_row(centers_mat, center);
    gsl_vector_sub(&ri.vector, cm);
    gsl_vector_set(mag2, center
		   ,gsl_pow_2(gsl_blas_dnrm2(&ri.vector)));
  }

  //digonals
  for (size_t alpha = 0; alpha < space; alpha++){
    double I_val = 0.0;
    for (size_t center = 0; center < centers; center++){
      I_val += gsl_vector_get(&mass.vector,center)	
	* (gsl_vector_get(mag2, center) - gsl_pow_2(gsl_matrix_get(centers_mat,center, alpha)));
    }
  gsl_matrix_set(I, alpha, alpha, I_val);
  }

  //off diags
  for (size_t row = 0; row < space; row++){
    for (size_t col = 0; col < row; col ++){
      double I_val = 0.0;
      for (size_t center = 0; center < centers; center++){
	I_val += -1.0 * gsl_vector_get(&mass.vector, center)
	  * gsl_matrix_get(centers_mat,center, row)
	  * gsl_matrix_get(centers_mat,center, col); //row and col reffer to cols in centers_mat
      }
      gsl_matrix_set(I, row, col, I_val);
      gsl_matrix_set(I, col, row, I_val);
    }
  }

  gsl_vector_free(cm);
  gsl_vector_free(mag2);
  gsl_matrix_free(centers_mat);
}

void get_angular_momentum(const gsl_vector * r, const gsl_vector * m, const gsl_vector * v, gsl_vector * L){
  const size_t space = 3; // gives the dimensionality of the space
  size_t dim = v->size;
  assert ((dim == r->size) && (dim == m->size));
  assert (L->size == space);
  assert (dim % space == 0);
  size_t N = dim / space;

  gsl_vector_set_zero(L);

  gsl_vector * cm = gsl_vector_alloc(space);
  center_of_mass(r, m, cm);

  gsl_vector_const_view mass = gsl_vector_const_subvector_with_stride(m, (size_t) 0, space, N);
  gsl_vector_view ri;
  gsl_vector_view vi;

  gsl_vector * li = gsl_vector_calloc(space);

  for (size_t i = 0; i < N; i++){
    // ri = gsl_vector_subvector(r, i * space, space);
    ri = gsl_vector_view_array(r->data + (i * space * r->stride),  space);
    // vi = gsl_vector_subvector(v, i * space, space);
    vi = gsl_vector_view_array(v->data + (i * space * v->stride),  space);
    gsl_vector_memcpy(li, &ri.vector);
    gsl_vector_sub(li, cm);
    vector_cross(li, &vi.vector);
    gsl_blas_daxpy(gsl_vector_get(&mass.vector,i), li, L);
  }

  gsl_vector_free(li);
  gsl_vector_free(cm);
}

// inverts the 3x3 matrix M in-place
// uses analytical inversion from https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3.C3.973_matrices
void invert_3x3(gsl_matrix * M){
  assert(M -> size1 == 3 && M -> size2 == 3);
  gsl_matrix * t = gsl_matrix_calloc((size_t) 3, (size_t) 3);
  gsl_matrix_memcpy(t, M);

  double det = det_3x3(M);

  gsl_vector_view r1, r2, r3, c1, c2, c3;
  r1 = gsl_matrix_row(M, (size_t) 0);
  r2 = gsl_matrix_row(M, (size_t) 1);
  r3 = gsl_matrix_row(M, (size_t) 2);

  c1 = gsl_matrix_column(t, (size_t) 0);
  c2 = gsl_matrix_column(t, (size_t) 1);
  c3 = gsl_matrix_column(t, (size_t) 2);

  //[ 2 x 3]
  //[ 3 x 1]
  //[ 1 x 2]

  gsl_vector_memcpy(&r1.vector, &c2.vector);
  vector_cross     (&r1.vector, &c3.vector);

  gsl_vector_memcpy(&r2.vector, &c3.vector);
  vector_cross     (&r2.vector, &c1.vector);

  gsl_vector_memcpy(&r3.vector, &c1.vector);
  vector_cross     (&r3.vector, &c2.vector);

  gsl_matrix_scale(M,1.0/det);

  gsl_matrix_free(t);
}


// returns the determinant of a 3x3 matrix
double det_3x3(gsl_matrix * M){
  assert(M -> size1 == 3 && M -> size2 == 3);

  double a11, a12, a13, a21, a22, a23, a31, a32, a33; 
  a11 = gsl_matrix_get(M,(size_t) 0, (size_t) 0);
  a12 = gsl_matrix_get(M,(size_t) 0, (size_t) 1);
  a13 = gsl_matrix_get(M,(size_t) 0, (size_t) 2);
  a21 = gsl_matrix_get(M,(size_t) 1, (size_t) 0);
  a22 = gsl_matrix_get(M,(size_t) 1, (size_t) 1);
  a23 = gsl_matrix_get(M,(size_t) 1, (size_t) 2);
  a31 = gsl_matrix_get(M,(size_t) 2, (size_t) 0);
  a32 = gsl_matrix_get(M,(size_t) 2, (size_t) 1);
  a33 = gsl_matrix_get(M,(size_t) 2, (size_t) 2);

  return a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - (a11*a23*a32 + a22*a13*a31 + a33*a12*a21);
}

void center_of_mass(const gsl_vector * r, const gsl_vector * m, gsl_vector * o){
  const size_t space = 3; // gives the dimensionality of the space
  size_t dim = r->size;
  assert (dim == m->size);
  assert (o->size == space);
  assert (dim % space == 0);
  size_t N = dim / space;

  gsl_vector_set_zero(o);

  gsl_vector_const_view mass = gsl_vector_const_subvector_with_stride(m, (size_t) 0, space, N);
  double m_total = 0.0;
  double mi;
  gsl_vector_view ri;

  
  for (size_t i = 0; i < N; i++){
    //ri = gsl_vector_subvector(r, i * space, space);
    ri = gsl_vector_view_array(r->data + (i * space * r->stride),  space);
    mi = gsl_vector_get(&mass.vector, i);
    gsl_blas_daxpy(mi,&ri.vector,o);
    m_total += mi;
  }
  
  gsl_vector_scale(o, 1/m_total);
}

//compute CM motion
void center_of_mass_velocity(const gsl_vector * v, const gsl_vector * m, gsl_vector * o){
  const size_t space = 3; // gives the dimensionality of the space
  size_t dim = v->size;
  assert (dim == m->size);
  assert (o->size == space);
  assert (dim % space == 0);
  size_t centers = dim / space;

  double mtotal = 0.0; 
  gsl_vector_set_zero(o);

  for (size_t center = 0; center < centers; center++){
    gsl_vector_const_view vi = gsl_vector_const_subvector(v, 3 * center, (size_t) 3);//velocity of each center
    double mi = gsl_vector_get(m, center * 3);
    mtotal += mi;
    gsl_blas_daxpy(mi, &vi.vector, o);
  }

  gsl_vector_scale(o, 1.0/mtotal);
}

// a <- a x b
void vector_cross(gsl_vector * a, gsl_vector * b){
  const size_t dim = 3;
  assert ((a->size == dim) && (b->size == dim));

  double a1, a2, a3, b1, b2, b3;
  double s1, s2, s3;

  a1 = gsl_vector_get(a, (size_t) 0);
  a2 = gsl_vector_get(a, (size_t) 1);
  a3 = gsl_vector_get(a, (size_t) 2);
  b1 = gsl_vector_get(b, (size_t) 0);
  b2 = gsl_vector_get(b, (size_t) 1);
  b3 = gsl_vector_get(b, (size_t) 2);

  s1 = (a2 * b3) - (a3 * b2);
  s2 = (a3 * b1) - (a1 * b3);
  s3 = (a1 * b2) - (a2 * b1);

  gsl_vector_set(a, (size_t) 0, s1);
  gsl_vector_set(a, (size_t) 1, s2);
  gsl_vector_set(a, (size_t) 2, s3);
}

//returns the euclidian distance between x1 and x2
double ext_gsl_vector_dist(gsl_vector * x1, gsl_vector * x2){
  size_t dim = x1->size;
  assert(dim == x2->size);
  
  double r;
  gsl_vector * t = gsl_vector_alloc(dim);

  gsl_vector_memcpy(t,x2);
  gsl_vector_sub(t,x1);
  r = gsl_blas_dnrm2(t);
  gsl_vector_free(t);
  return r;
}
