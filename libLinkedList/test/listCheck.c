#include "listCheck.h"
#include "../src/linkedList.h"

struct pointerArray{
  size_t nElems;
  size_t size;
  void ** a;
};

static struct pointerArray dataStore = { .nElems = 0, .size = 0, .a = NULL};

static void nullCrash(void);
static void crashOnFail(llist * L);

static size_t store(void * p);
static size_t lookup(void * p);
static void clear(void);

// Initiates a segfault if the properties of a list are violated
// see 2014.09.09:184
void listCheck(llist * L){
  llist * t;
  clear();
  //forward walk
  for (t=L; t != NULL; t=t->succ){
    crashOnFail(t);
  }

  //reverse walk; start at pred so only touch each elem once
  for (t=L->pred; t != NULL; t=t->pred){
    crashOnFail(t);
  }

}

static void crashOnFail(llist * L){
  if (!L)
    return;

  if (lookup(L)) //found a dupe!
    nullCrash();
  else
    store(L);

  if ((L->succ != NULL) && (L != L->succ->pred))
    nullCrash();

  if ((L->pred != NULL) && (L != L->pred->succ))
    nullCrash();
}

// Causes a segfault and the writing of a core dump
static void nullCrash(void){
  *((int *) 0) = 0;
}

static void clear(void){
  dataStore.nElems = 0;
}

static size_t store(void * p){
  const size_t increase = 100;
  if (dataStore.nElems == dataStore.size){
    void ** t = realloc(dataStore.a
			, (dataStore.size + increase) * sizeof(void *));
    if (t){
      dataStore.a = t;
      dataStore.size += increase;
    }
    else
      return 0;
  }

  dataStore.a[dataStore.nElems] = p;
  return dataStore.nElems++;
}

static size_t lookup(void * p){
  for (size_t i = 0; i < dataStore.nElems; i++){
    if (p == dataStore.a[i])
      return 1;
  }
  return 0;
}
