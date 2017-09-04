#include "kinds.h"
#include <stdio.h>

#define BUF 1024

/*
  The memory returned is static and should not be freed
*/
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


geodesicKind c2k(char c){
  switch (c){
  case 'D': /* exclicit fallthrough */
  case 'd':
    return direct;
    break;
  case 'R': /* exclicit fallthrough */
  case 'r':
    return roaming;
    break;
  case 'A': /* exclicit fallthrough */
  case 'a':
    return radical;
    break;
  case 'N': /* exclicit fallthrough */
  case 'n':
  default:
    return nothing;
    break;
  }
}
