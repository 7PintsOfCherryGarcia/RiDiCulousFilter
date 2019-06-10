#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "khash.h"
#include "kseq.h"
#include "filters.h"


int usage(char **argv) {
  fprintf(stderr,"Usage: %s %s -l <int> -u <int> -f FILE\n", argv[0], argv[1]);
  return -1;
}

int mainCountFilter(int argc, char **argv) {
  int elem;
  opts opt;
  opt.file = NULL;
  while (( elem = getopt(argc, argv, ":l:u:f:") ) >= 0) {
    switch(elem) {
    case 'l':
      opt.lower = atoi(optarg);
      break;
    case 'u':
      opt.upper = atoi(optarg);
      break;
    case 'f':
      opt.file = optarg; break;
    }
  }

  if (!opt.lower || !opt.upper || !opt.file) {
    return usage(argv);
  }
  printf("%d, %d, %s\n",opt.lower,opt.upper, opt.file);


  //Read arguments, filter type, sequence files
  FILE *fp;
  // Read from stdin
  if(strcmp(opt.file,"-") == 0) {
    opt.file = "/dev/stdin";
  }
  fp = fopen(opt.file, "r");
  if (!fp) {
    fprintf(stderr,"Can't open input file %s.\n", opt.file);
    return(usage(argv));
  }

  int numKmers = 0;
  kmerCount *counts = malloc(10000*sizeof(kmerCount));
  if (!counts) {
    fprintf(stderr, "Error insufficient memeory.\n");
    return -1;
  }
  counts = readKmerCounts(fp, counts, &numKmers);
  if(!counts) {
    fprintf(stderr, "Error loading kmer list\n");
    return -1;
  }
  fprintf(stderr,"%d kmers\n", numKmers);
  fclose(fp);

  //Load filter into hash table
  //Create and initialize hash table
  //TODO free khash varables
  khash_t(kmer) *h;
  h = kh_init(kmer);
  khint_t k;
  k = kh_end(h);
  int ret = loadKmerHash(h, k, counts, 10);
  return ret;
  //Loop over reads apply filter and write to stdout
}

kmerCount *readKmerCounts(FILE *fp, kmerCount *counts, int *numKmers) {
  printf("Reading\n");
  char *line = NULL;
  size_t n = 0;
  char *kmer;
  int count;
  char *err;
  int numChars;

  while((numChars = getline(&line, &n, fp)) != -1) {
    line[numChars - 1] = '\0';
    kmer = strtok(line, "\t");
    count = strtol(strtok(NULL,"\t"), &err, 10);
    if(*err != 0) {
      fprintf(stderr,"Error in line %d\n",*numKmers + 1);
      return NULL;
    }
    strcpy(counts[*numKmers].kmer,kmer);
    counts[*numKmers].count = count;

    *numKmers += 1;
    if(*numKmers%1000 == 0) {
      counts = (kmerCount *) realloc(counts,(*numKmers + 10000)*sizeof(kmerCount));
    }
  }

  free(line);
  free(counts);
  fprintf(stderr,"%d kmers\n", *numKmers);
  return counts;
}

int loadKmerHash(khash_t(kmer) *h, khint_t k, kmerCount *counts, int numKmers) {
  
}
