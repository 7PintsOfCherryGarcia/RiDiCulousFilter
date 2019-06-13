#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>


int main_CountFilter(int, char **);
int mainGCFilter(int argc, char **argv) {
  return 0;
}

int main_usage(char **argv) {
  fprintf(stderr,"Usage: %s [command] [options] -f FILE\n", argv[0]);
  return -1;
}

int main(int argc, char *argv[]) {

  if(argc <= 2) {
    return main_usage(argv);
  }

  if(strcmp(argv[1],"count") == 0) {
    fprintf(stderr,"Filtering by count\n");
    return main_CountFilter(argc - 1, argv + 1);
  }
  else if(strcmp(argv[1],"GC") == 0) {
    fprintf(stderr,"Filtering by GC content\n");
    return mainGCFilter(argc - 1, argv + 1);
  }
  else {
    fprintf(stderr,"Command: \"%s\" not recognized.\n",argv[1]);
    main_usage(argv);
  }
  return(-1);
}
