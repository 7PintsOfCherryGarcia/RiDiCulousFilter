#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "../khash.h"
#include "../kseq.h"


typedef struct {
  char *kmer;
  int count;
} kmerCount;

int readKmerCounts(FILE *, kmerCount *);

int main(int argc, char *argv[]) {
  //Read arguments, filter type, sequence files
  FILE *fp;
  // Read from stdin
  fp = fopen(argv[2], "r");
  if (!fp) {
    printf("Can't open input file.\n");
    exit(1);
    //usage();
  }

  kmerCount *counts = malloc(10000*sizeof(kmerCount));
  int tmp = readKmerCounts(fp, counts);
  printf("%d kmers parsed\n",tmp);
  fclose(fp);
  free(counts);
  //Load filter into hash table

  //Loop over reads apply filter and write to stdout


  return(0);
}
