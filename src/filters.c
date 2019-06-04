#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>

#include "../khash.h"
#include "../kseq.h"
#include "filters.h"




int readKmerCounts(FILE *fp, kmerCount *counts) {
  printf("Reading\n");
  char *line = NULL;
  size_t n = 0;
  int kmerCounter = 0;
  char *kmer;
  char *count;
  while(getline(&line, &n, fp) != -1) {
    kmer = strtok(line, "\t");
    printf("%s\n",kmer);
    count = strtok(NULL,"\t");
    printf("\t%s",count);
    kmerCounter += 1;
  }
  free(line);
  return(kmerCounter);
}
