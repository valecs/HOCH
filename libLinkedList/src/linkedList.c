#include "linkedList.h"
#include <stdio.h>

#include <assert.h>

/*
TODO:
* eliminate asserts
* rewrite with for- instead of while- loops
* implement stressful test framework
 */

llist * ll_new(void * elem){
  llist * head = malloc(sizeof(llist));
  head -> pred = NULL;
  head -> succ = NULL;
  head -> e = elem;
  return head;
}

/* void ll_destroy(llist * ) */

//you can always append(NULL,elem) to get started;
llist * append(llist * curr, void * elem){
  llist * butt = ll_new(elem);
  if ( !curr ) {return butt;}//empty list
  //get last element
  while (curr -> succ != NULL) {curr = curr -> succ;}
  curr -> succ = butt;
  butt -> pred = curr;
  return butt;
}

/* count both ways from starting point then add one for initial */
unsigned int length(llist * curr){
  unsigned int down = 0, up = 0;

  if ( !curr ){return 0;} //empty list

  llist * init = curr;
  while (curr -> succ != NULL) {curr = curr -> succ; up++;}
  curr=init;
  while (curr -> pred != NULL) {curr = curr -> pred; down++;}
  return 1 + down + up;
}

llist * last(llist * curr){
  if ( !curr ) { return NULL; }//empty list
  while ((curr -> succ) != NULL)
    {curr = curr -> succ;}
  return curr;
}

llist * first(llist * curr){
  if ( !curr ) { return NULL; }//empty list
  while (curr -> pred != NULL) {curr = curr -> pred;}
  return curr;
}

//doesn't dealocate the element
//returns the list without it's head, dealocating it
llist * tail(llist * curr){
  if ( !curr ) {return NULL;} //empty list
  
  curr = first(curr);
  if (curr->succ){
    curr=curr->succ;
    free(curr->pred);
    curr->pred=NULL;
    return curr;
  } else{
    free(curr);
    return NULL;
  }

}

//reverse a list!
llist * reverse(llist * path){
  llist * priorFirst = first(path);
  llist * priorLast = last(path);
  llist * swap;
  for (llist * L = first(path); L != NULL; L = L->pred){ // L=L->pred :: pseudo reverse walk
    swap = L->pred;
    L->pred = L->succ;
    L->succ = swap;
  }

  assert(last(path) == priorFirst);
  assert(first(path) == priorLast);
  return first(path);
}

//returns a union b with simple loop checking
llist * join(llist * a, llist * b){
  if (!a && !b) {return NULL;} //the union of empty lists is empty
  assert(first(a)!=first(b) && last(a) != last(b));
  return join_unsafe(a,b);
}

//returns a union b
//no checking for loops
llist * join_unsafe(llist * a, llist * b){
  if (!a && !b) {
    return NULL;//the union of empty lists is empty
  } else if (a && !b) {
    return a;
  } else if (b && !a) {
    return b;
  } else{
    llist * a_end = last(a);
    a_end->succ=first(b);
    first(b)->pred=a_end;
    return a;
  }
}

// returns pointer to ith element; indexes start at 0.
// no bounds checking
llist * lindex_unsafe(llist * curr, unsigned int ind){
  if ( !curr ) { return NULL; }//empty list
  curr = first(curr);
  for (unsigned int i = 0; i < ind ; i++){
    curr = curr->succ;
  }
  return curr;
}

// list indexing with buounds checking
llist * lindex(llist * curr, unsigned int ind){
  if ( !curr ) { return NULL; }//empty list
  assert(ind<length(curr));
  return lindex_unsafe(curr, ind);
}

/* void * fold1(llist * curr, (void *) (*f)(void *)){ */
/*   curr = first(curr); */
/* } */

llist * map(llist * l, mapFunction_t f){
  void * e;

  if ( !l ){return NULL;} //return empty list if list is empty

  llist * result = NULL; //the empty list

  l = first(l);
  
  while ( l->succ){
    e = (*f)(l->e);
    if (e){
      result = append(result,e);
    }
    l=l->succ;
  }

  //catching the last element
  e = (*f)(l->e);
  if (e){
    result = append(result,e);
  }
  
  return result;
}
