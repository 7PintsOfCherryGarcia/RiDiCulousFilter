#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "khash.h"
#include "kseq.h"
#include "filters.h"

static int true = 1;

int usage(char **argv) {
  fprintf(stderr,"Usage:\nRiDiCulous %s -l <int> -u <int> -c KMER_FILE -f SEQ_FILE\n", argv[0]);
  return -1;
}

int main_CountFilter(int argc, char **argv) {
  int elem;
  opts opt;
  opt.kmerfile = NULL;
  opt.seqfile = NULL;
  opt.seqFP = NULL;
  opt.minfraction = 0.80;
  opt.kmerlen = 31;
  while (( elem = getopt(argc, argv, "l:u:c:f:m:k:") ) >= 0) {
    switch(elem) {
    case 'l':
      opt.lower = atoi(optarg);
      break;
    case 'u':
      opt.upper = atoi(optarg);
      break;
    case 'c':
      opt.kmerfile = optarg;
      break;
    case 'f':
      opt.seqfile = optarg;
      break;
    case 'm':
      opt.minfraction = atof(optarg);
      break;
    case 'k':
      opt.kmerlen = atoi(optarg);
    }
  }

  if (!opt.lower || !opt.upper || !opt.kmerfile || !opt.seqfile) {
    return usage(argv);
  }

  if(opt.minfraction <= 0 || opt.minfraction >= 1) {
    opt.minfraction = 0.80;
  }
  fprintf(stderr,"%d, %d, %d, %s, %s\n",
         opt.lower,opt.upper, opt.kmerlen, opt.kmerfile, opt.seqfile);
  fprintf(stderr,"%f\n",opt.minfraction);

  //Read arguments, filter type, sequence files
  FILE *fp;
  int stdinFlag = 0;
  // Read from stdin
  if(strcmp(opt.kmerfile,"-") == 0) {
    opt.kmerfile = "/dev/stdin";
    stdinFlag = 1;
  }
  fp = fopen(opt.kmerfile, "r");
  if (!fp) {
    fprintf(stderr,"Can't open input file %s.\n", opt.kmerfile);
    return(usage(argv));
  }

  if(strcmp(opt.seqfile,"-") == 0) {
    opt.seqfile = "/dev/stdin";
    if(stdinFlag == 1) {
      fprintf(stderr,"Error: Only one input file can be set to stdin.");
      usage(argv);
      return -1;
    }
  }
  opt.seqFP = gzopen(opt.seqfile,"r");
  if(!opt.seqFP) {
    fprintf(stderr,"Error: Can't open sequence file.");
    return -1;
  }



  //Khash
  khash_t(kmer) *h;
  h = kh_init(kmer);
  khint_t k;
  k = kh_end(h);

  fprintf(stderr,"Loading kmers\n");
  int kmerCount = readKmerCounts(fp, opt.kmerlen, h, &k);
  if(kmerCount <= 0) {
    fprintf(stderr, "Error loading kmer list\n");
    return -1;
  }
  fprintf(stderr,"%d kmers\n", kmerCount);
  fclose(fp);


  fprintf(stderr,"Loaded %d kmer counts.\n",kmerCount);

  //Loop over reads apply filter and write to stdout
  fprintf(stderr,"Filtering reads.\n");
  unsigned long numReads = filterReads(h, k, opt);
  fprintf(stderr,"%ld reads processed\n",numReads);
  hash_destroy(h, k);
  gzclose(opt.seqFP);
  return 0;
}


int readKmerCounts(FILE *fp,
                   int kmerlen,
                   khash_t(kmer) *h,
                   khint_t *k) {

  int kmerCount = 0;
  char *line = NULL;
  size_t n = 0;
  char *kmer;
  char *countStr;
  int count;
  char *err;
  int numChars;
  int lineNum = 1;
  int absent = 0;

  while((numChars = getline(&line, &n, fp)) != -1) {
    line[numChars - 1] = '\0';
    kmer = strtok(line, "\t");
    if(!kmer) {
      fprintf(stderr,
              "WARNING: In line %d\nCould not read kmer sequence. Skipping\n",
              lineNum);
      lineNum += 1;
      continue;
    }
    else if(strlen(kmer) != kmerlen) {
      fprintf(stderr,
              "WARNING in line %d\nkmer length of different size. Skipping\n",
              lineNum);
      lineNum += 1;
      continue;
    }
    countStr = strtok(NULL,"\t");
    if(!countStr) {
      fprintf(stderr,
              "WARNING: In line %d\nkmer %s without a count. Skipping\n",
              lineNum, kmer);
      lineNum += 1;
      continue;
    }
    count = strtol(countStr, &err, 10);
    if(*err != 0) {
      fprintf(stderr,
              "WARNING in line %d\ncould not read kmer count. Skipping\n",
              lineNum);
      lineNum += 1;
      continue;
    }


    //Load into hash
    *k = kh_put(kmer, h, kmer,&absent);
    kh_key(h, *k) = strdup(kmer);
    kh_value(h,*k) = count;


    kmerCount += 1;
    lineNum += 1;

  }

  free(line);
  return kmerCount;
}


//unsigned int loadKmerHash(khash_t(kmer) *h,
                          //                          khint_t k,
                          //                          kmerCount *counts,
                          //                          int numKmers) {
  //  int absent = 0;
  //  for(int i = 0; i < numKmers; i++) {
    //    k = kh_put(kmer, h, counts[i].kmer,&absent);
    //    kh_key(h, k) = strdup(counts[i].kmer);
    //   kh_value(h,k) = counts[i].count;
    //  }
  //  return kh_size(h);
  //}


unsigned long filterReads(khash_t(kmer) *h,
                          khint_t k,
                          opts opt) {
  kseq_t *seq;
  int l;
  int n = 0;
  unsigned long int totalSeq = 0;
  char *kmer = malloc((opt.kmerlen + 1)*sizeof(char));
  kmer[opt.kmerlen] = '\0';

  seq = kseq_init(opt.seqFP);
  while ((l = kseq_read(seq)) >= 0) {
    if (l == 0) continue;
    totalSeq += l;
    if (l <= 0) {
      fprintf(stderr,"Error in sequence file.\n");
      return -1;
    }
    if(query_read(h, k, opt, seq->seq.s, l, kmer)) {
      printf("@%s\n%s\n+\n%s\n",seq->name.s,seq->seq.s,seq->qual.s);
    }
    else {
      //fprintf(stderr,"@%s\n%s\n+\n%s\n",seq->name.s,seq->seq.s,seq->qual.s);
    }

    n++;
  }
  kseq_destroy(seq);
  free(kmer);
  return n;
}


int *query_read(khash_t(kmer) *h,
                khint_t k,
                opts opt,
                char *seq,
                int l,
                char *kmer) {
  int kmernum = 0;
  double kmerFraction;
  for(int i = 0; i < l - opt.kmerlen + 1; i++) {
    strncpy(kmer, &seq[i], opt.kmerlen);
    k = kh_get(kmer,h,kmer);
    if (k == kh_end(h)) {
      continue;
    }
    else if((kh_value(h,k) >= opt.lower) && (kh_value(h,k) <= opt.upper)) {
      kmernum += 1;
    }
  }
  kmerFraction = (double) kmernum / (double) (l - opt.kmerlen + 1);
  if (kmerFraction >= opt.minfraction) {
    return &true;
  }
  else {
    return NULL;
  }
}


void hash_print(khash_t(kmer) *h, khint_t k) {
  //Loop over hash integers and print key value pair
  for (k = 0; k < kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      printf("%s\t%d\n",kh_key(h,k),kh_val(h,k));
    }
  }
}


void hash_destroy(khash_t(kmer) *h, khint_t k) {
  for (k = 0; k < kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      free((char*)kh_key(h, k));
    }
  }
  kh_destroy(kmer, h);
}
