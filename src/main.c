#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>


int main_count(int, char **);
int main_gc(int argc, char **argv) {
  return 0;
}
int main_kmer(int argc, char **argv){
  return 0;
}

static int main_usage() {
  fprintf(stderr,"Usage: RiDiCulous [command] [options] -f FILE\n");
  fprintf(stderr,"Commands:\n");
  fprintf(stderr,"\tcount -  Filter reads based on kmer counts.\n\n");
  fprintf(stderr,"\t         You must suply a kmer count table.\n\n");
  fprintf(stderr,"\tkmer  -  Filter reads based on kmer presence.\n\n");
  fprintf(stderr,"\t         You must suply a kmer table.\n\n");
  fprintf(stderr,"\tgc    -  Filter reads based on GC content.\n\n");
  fprintf(stderr,"Help information for these commands is diplayed on\n");
  fprintf(stderr,"erroneous execution or when the -h flag is provided.\n\n");
  return -1;
}

int main(int argc, char *argv[]) {

  if(argc <= 1) {
    return main_usage();
  }

  if(strcmp(argv[1],"count") == 0) {
    fprintf(stderr,"Command: count\n");
    return main_count(argc - 1, argv + 1);
  }
  else if(strcmp(argv[1],"gc") == 0) {
    fprintf(stderr,"command: gc\n");
    return main_gc(argc - 1, argv + 1);
  }
  else if(strcmp(argv[1],"kmer") == 0) {
    fprintf(stderr,"command: kmer\n");
    return main_kmer(argc - 1, argv + 1);
  }
  else {
    fprintf(stderr,"Command: \"%s\" not recognized.\n",argv[1]);
    main_usage();
  }
  return(-1);
}
