#include "hcohHelpers.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <assert.h>

//takes a total phase space vector and returns them for each product molecule
void get_product_vectors(gsl_vector * g, gsl_vector * m
			 ,gsl_vector * g_HH, gsl_vector * m_HH
			 ,gsl_vector * g_CO, gsl_vector * m_CO){

  const size_t dim    = 12;
  const size_t dim_HH =  6;
  const size_t dim_CO =  6;
  const size_t space  =  3;
  
  assert(g->size    == 2 * dim);    assert(m->size    == dim);
  assert(g_HH->size == 2 * dim_HH); assert(m_HH->size == dim_HH);
  assert(g_CO->size == 2 * dim_CO); assert(m_CO->size == dim_CO);
  
  /*
  N.B.: all vectors in the simulation are in the order H1, O, C, H2

  on return HH: H1, H2
            CO: C,  O

  */
  //chunks for moving things around
  gsl_vector_view rxn, prod;

  //the plan is to set rxn & prod then memcpy() for a long list of combinations.

  //masses
  rxn  = gsl_vector_subvector(m,    0 * space, space); //H1
  prod = gsl_vector_subvector(m_HH, 0 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(m,    3 * space, space); //H2
  prod = gsl_vector_subvector(m_HH, 1 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(m,    1 * space, space); //O
  prod = gsl_vector_subvector(m_CO, 1 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(m,    2 * space, space); //C
  prod = gsl_vector_subvector(m_CO, 0 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  //positions
  rxn  = gsl_vector_subvector(g,    0 * space, space); //H1
  prod = gsl_vector_subvector(g_HH, 0 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(g,    3 * space, space); //H2
  prod = gsl_vector_subvector(g_HH, 1 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(g,    1 * space, space); //O
  prod = gsl_vector_subvector(g_CO, 1 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(g,    2 * space, space); //C
  prod = gsl_vector_subvector(g_CO, 0 * space, space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  //velocities
  rxn  = gsl_vector_subvector(g,    dim    + (0 * space), space); //H1
  prod = gsl_vector_subvector(g_HH, dim_HH + (0 * space), space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(g,    dim    + (3 * space), space); //H2
  prod = gsl_vector_subvector(g_HH, dim_HH + (1 * space), space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(g,    dim    + (1 * space), space); //O
  prod = gsl_vector_subvector(g_CO, dim_CO + (1 * space), space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

  rxn  = gsl_vector_subvector(g,    dim    + (2 * space), space); //C
  prod = gsl_vector_subvector(g_CO, dim_CO + (0 * space), space);
  gsl_vector_memcpy(&prod.vector, &rxn.vector); //dest,src

}

// returns CO <-> H2 CM distance in a0
//should not touch contents of x or mass
double dist_H2_CO(gsl_vector * x, gsl_vector * mass){
  static const size_t dim = 3;
  
  static gsl_vector * rCO = NULL;
  static gsl_vector * d   = NULL;
  
  if (rCO == NULL){rCO = gsl_vector_alloc(dim);}
  if (d   == NULL){d   = gsl_vector_alloc(dim);}

  gsl_vector_view rH1 = gsl_vector_subvector(x, 0 * dim, dim);
  gsl_vector_view rO  = gsl_vector_subvector(x, 1 * dim, dim);
  gsl_vector_view rC  = gsl_vector_subvector(x, 2 * dim, dim);
  gsl_vector_view rH2 = gsl_vector_subvector(x, 3 * dim, dim);

  double mO = gsl_vector_get(mass, 1 * dim);
  double mC = gsl_vector_get(mass, 2 * dim);

  //use d as a temp var:
  gsl_vector_memcpy(d, &rC.vector);
  gsl_vector_memcpy(rCO, &rO.vector);
  gsl_vector_scale( rCO, mO);
  gsl_blas_daxpy(mC, d, rCO);//y<-ax+y
  gsl_vector_scale(rCO, (1.0/(mC+mO)));

  gsl_vector_memcpy(d, &rH1.vector);
  gsl_vector_add(   d, &rH2.vector);
  gsl_vector_scale( d, 0.5);

  gsl_vector_sub(d, rCO);  

  return gsl_blas_dnrm2(d);
}


void get_hoch_center_distances(gsl_vector * x, gsl_vector * r){
  const size_t space = 3;
  assert (x -> size % space == 0);
  size_t N = x->size / 3; 
  assert (r -> size == N*(N-1)/2);

  size_t idx = 0;

  gsl_vector_view a,b;
  gsl_vector * t = gsl_vector_alloc(space);

  for (size_t i = 0; i < N; i++){
    a = gsl_vector_subvector(x, i * space, space);
    for(size_t j = 0; j < i; j++){
      gsl_vector_memcpy(t,&a.vector);
      b = gsl_vector_subvector(x, j * space, space);
      gsl_vector_sub(t, &b.vector);
      gsl_vector_set(r, idx, gsl_blas_dnrm2(t));
      idx++; //idx = i + j + [i(i-3)/2], but I presume this is faster
    }
  }

  gsl_vector_free(t);
}
