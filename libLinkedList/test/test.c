#include "../src/linkedList.h"
#include "listCheck.h"

#include <stdio.h>

static llist * construct (int n);

int main(void){
  const int len = 5;
  llist * L = construct(len);
  printf("Constructed list of length %d==%d\n", length(L), len);

  listCheck(L);
  printf("listCheck passed!\n");

  //now we break!
  //L = first(L);
  //L->succ->succ->pred=L;

  listCheck(last(L));
  printf("listCheck passed!\n");

  return 0;
}

static llist * construct (int n){
  int e = 0;
  llist * L = NULL;
  for (int i = 0; i<n; i++){
    L=append(L,&e);
  }
  return L;
}
