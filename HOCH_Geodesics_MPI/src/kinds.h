#ifndef KINDS_H
#define KINDS_H

typedef enum {
  nothing,
  direct,
  roaming,
  radical
} geodesicKind;

const char * k2s(geodesicKind k);
char k2c(geodesicKind k);
geodesicKind c2k(char c);

#endif
