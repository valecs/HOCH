#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <stdlib.h>

typedef struct llist llist;
struct llist { 
  llist * pred;
  llist * succ;
  void * e;
};

typedef void * (*mapFunction_t)(void *);

llist * ll_new(void * elem);
llist * append(llist * curr, void * elem);
unsigned int length(llist * curr);
llist * last(llist * curr);
llist * first(llist * curr);
llist * tail(llist * curr);
llist * reverse(llist * path);
llist * join(llist * a, llist * b);
llist * join_unsafe(llist * a, llist * b);
llist * lindex(llist * curr, unsigned int ind);
llist * lindex_unsafe(llist * curr, unsigned int ind);
llist * map(llist * l, mapFunction_t f);

#endif
