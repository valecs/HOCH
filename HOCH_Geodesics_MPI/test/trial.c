#include <stdio.h>
#include <stdlib.h>

#define BUF 7

void msg(const char * s);
void getTaskFromFile(const char * fname);

int manifestID = -1;

int main(int argc, char ** argv){
  if (argc != 3){
    msg("Require 2 arguments.");
    exit(1);
  }

  manifestID = atoi(argv[2]);
  fprintf(stderr, "Got: %d\n", manifestID);

  getTaskFromFile(argv[1]);

  printf("manifestID = %05d\n", manifestID);

  return 0;
}

void msg(const char * s){
  fprintf(stderr, "%s\n", s);
}

void getTaskFromFile(const char * fname){
  FILE * f;
  f = fopen(fname, "r");
  if (!f){
    msg(fname);
    msg("FAILURE_READ_OPEN");
    exit(1);
  }
  
  char in[BUF];
  char out[BUF];
  int j = 0;
  int lines = manifestID;
  while (0 != fread(in, sizeof(char), BUF, f) && lines > 0){
    for(char * i = &in[0]; i < &in[BUF]; i++){
      if (1 == lines){
	out[j++] = (*i);
      }

      if (j >= BUF){
	msg("length error");
	return;
      }

      if ('\n' == (*i) && !--lines){
	  out[j-1] = '\0';
	  manifestID = atoi(out);
	  break;
      }
    }
  }

  fclose(f);
}
