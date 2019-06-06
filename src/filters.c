#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "../khash.h"
#include "../kseq.h"
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

  kmerCount *counts = malloc(10000*sizeof(kmerCount));
  counts = readKmerCounts(fp, counts);

  fclose(fp);
  free(counts);
  //Load filter into hash table

  //Loop over reads apply filter and write to stdout
  return(0);
}

kmerCount *readKmerCounts(FILE *fp, kmerCount *counts) {
  printf("Reading\n");
  char *line = NULL;
  size_t n = 0;
  int kmerCounter = 0;
  char *kmer;
  int count;
  char *err;
  int numChars;
  while((numChars = getline(&line, &n, fp)) != -1) {
    line[numChars - 1] = '\0';
    kmer = strtok(line, "\t");
    count = strtol(strtok(NULL,"\t"), &err, 10);
    if(*err != 0) {
      fprintf(stderr,"Error in line %d\n",kmerCounter + 1);
      return NULL;
    }
    strcpy(counts[kmerCounter].kmer,kmer);
    counts[kmerCounter].count = count;
    //printf("%s\t",counts[kmerCounter].kmer);
    //printf("%d\n",counts[kmerCounter].count);

    kmerCounter += 1;
    if(kmerCounter%1000 == 0) {
      counts = (kmerCount *) realloc(counts,(kmerCounter + 10000)*sizeof(kmerCount));
    }
  }
  free(line);
  printf("EEEE\n");
  printf("%d kmers\n",kmerCounter);
  return(counts);
}
