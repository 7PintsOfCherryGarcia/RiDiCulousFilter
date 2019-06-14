#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "khash.h"
#include "kseq.h"
#include "main_count.h"

static int true = 1;


int count_usage() {
  fprintf(stderr,
          "Usage: RiDiCulous [options] -c KMER_FILE -f SEQ_FILE\n");
  fprintf(stderr,"Options:\n\n");
  fprintf(stderr,"\t-c\tkmer count file. Tab separated text file with\n");
  fprintf(stderr,"\t  \tkmer sequences as first column and counts in\n");
  fprintf(stderr,"\t  \tsecond column.\n\n");
  fprintf(stderr,"\t-f\tSequence file. fasta/q file with sequences or\n");
  fprintf(stderr,"\t  \tsequencing reads. Can be gziped. provide \"-\" if\n");
  fprintf(stderr,"\t  \treading from stdin. If using sequencing data eg.\n");
  fprintf(stderr,"\t  \tillumina reads, provide the -C flag to use canonical\n");
  fprintf(stderr,"\t  \tform.\n\n");
  fprintf(stderr,"\t-k\tkmer length. Different size kmers will be ignored.\n\n");
  fprintf(stderr,"\t-l\tminimum count threshold. Only kmers with counts\n");
  fprintf(stderr,"\t  \tgreater or equal to this value will be consifered.\n\n");
  fprintf(stderr,"\t-m\tMinimum coverage thershold. Sequences must be\n");
  fprintf(stderr,"\t  \tcovered by at least this fraction of count passing\n");
  fprintf(stderr,"\t  \tkmers. Otherwise, sequence is filtered out\n\n");
  fprintf(stderr,"\t-u\tmaximum count threshold. Only kmers with counts\n");
  fprintf(stderr,"\t  \tless than or equal to this value will be consifered.\n");
  fprintf(stderr,"\n");
  return -1;
}


void count_readOpt(int argc, char **argv, opts* opt) {
  int elem;
  while (( elem = getopt(argc, argv, "c:f:hk:l:m:u:") ) >= 0) {
    switch(elem) {
    case 'l':
      opt->lower = atoi(optarg);
      if(!opt->lower) {
        fprintf(stderr,"\tERROR: Option \"-l %s\" not valid\n\n",optarg);
        exit(count_usage());
      }
      break;
    case 'u':
      opt->upper = atoi(optarg);
      if(!opt->upper) {
        fprintf(stderr,"\tERROR: Option \"-u %s\" not valid\n\n",optarg);
        exit(count_usage());
      }
      break;
    case 'c':
      opt->kmerfile = optarg;
      break;
    case 'f':
      opt->seqfile = optarg;
      break;
    case 'h':
      exit(count_usage());
    case 'm':
      opt->minfraction = atof(optarg);
      if(!opt->minfraction) {
        fprintf(stderr,"\tERROR: Option \"-m %s\" not valid\n\n",optarg);
        exit(count_usage());
      }
      break;
    case 'k':
      opt->kmerlen = atoi(optarg);
    }
  }

  if (!opt->kmerfile || !opt->seqfile) {
    exit(count_usage());
  }

  if(opt->minfraction <= 0 || opt->minfraction > 1) {
    fprintf(stderr,"\tWARNING: minimum fraction not in \(0,1] range.\n");
    fprintf(stderr,"\t         Setting to default -m 0.8\n\n");
    opt->minfraction = 0.80;
  }

  count_printOpt(*opt);
}


void count_printOpt(opts opt) {
  fprintf(stderr,"\tkmer count file: %s\n",opt.kmerfile);
  fprintf(stderr,"\tsequence file: %s\n",opt.seqfile);
  fprintf(stderr,"\tlower count bound: %d\n",opt.lower);
  fprintf(stderr,"\tupper count bound: %d\n",opt.upper);
  fprintf(stderr,"\tminimum read coverage fraction: %f\n\n",opt.minfraction);
}


int main_count(int argc, char **argv) {
  //Set options
  opts opt;
  opt.kmerfile = NULL;
  opt.seqfile = NULL;
  opt.seqFP = NULL;
  opt.minfraction = 0.80;
  opt.kmerlen = 31;

  //Read options
  count_readOpt(argc, argv, &opt);

  //Open FILE options
  //Open kmer count file
  FILE *fp;
  int stdinFlag = 0;
  // Read from stdin
  if(strcmp(opt.kmerfile,"-") == 0) {
    opt.kmerfile = "/dev/stdin";
    stdinFlag = 1;
  }
  fp = fopen(opt.kmerfile, "r");
  if (!fp) {
    fprintf(stderr,"\tCan't open input file %s.\n", opt.kmerfile);
    return(count_usage());
  }

  //Open sequence file
  if(strcmp(opt.seqfile,"-") == 0) {
    opt.seqfile = "/dev/stdin";
    if(stdinFlag == 1) {
      fprintf(stderr,"\tERROR: -f - and  -c - are incompatibe.\n");
      fprintf(stderr,"\tOnly one input file can be read from stdin.\n\n");
      return count_usage();
    }
  }
  opt.seqFP = gzopen(opt.seqfile,"r");
  if(!opt.seqFP) {
    fprintf(stderr,"\tERROR: Can't open sequence file.\n\n");
    return count_usage();
  }



  //Declare hash objects
  khash_t(kmer) *h;
  h = kh_init(kmer);
  khint_t k;
  k = kh_end(h);

  //Load kmer counts into hash table
  fprintf(stderr,"\tLoading kmer count table\n");
  int kmerCount = 0;
  int rkCount = 0;
  count_readKmers(fp, opt, h, &k, &kmerCount, &rkCount);
  if(kmerCount <= 0) {
    fprintf(stderr, "\tERROR: kmer count table could not be parsed.\n\n");
    return count_usage();
  }
  fclose(fp);
  fprintf(stderr,"\tLoaded %d kmers.\n",kmerCount);
  fprintf(stderr,"\tSkipped %d kmers.\n",rkCount);


  //Loop over reads apply filter and write to stdout
  fprintf(stderr,"\tFiltering reads.\n");
  unsigned long numReads = 0;
  unsigned long numFiltered = 0;
  count_filterReads(h, k, opt, &numReads, &numFiltered);
  fprintf(stderr,"\t%ld reads processed\n",numReads);
  fprintf(stderr,"\t%ld reads filtered\n",numFiltered);
  fprintf(stderr,"\t%ld reads kept\n",numReads-numFiltered);
  hash_destroy(h, k);
  gzclose(opt.seqFP);
  return 0;
}


void count_readKmers(FILE *fp,
                    opts opt,
                    khash_t(kmer) *h,
                    khint_t *k,
                    int *kmerCount,
                    int *rkCount) {

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
              "\tWARNING: In line %d\nCould not read kmer sequence. Skipping\n",
              lineNum);
      lineNum += 1;
      continue;
    }
    else if(strlen(kmer) != opt.kmerlen) {
      fprintf(stderr,
              "\tWARNING in line %d\nkmer length of different size. Skipping\n",
              lineNum);
      lineNum += 1;
      continue;
    }
    countStr = strtok(NULL,"\t");
    if(!countStr) {
      fprintf(stderr,
              "\tWARNING: In line %d\nkmer %s without a count. Skipping\n",
              lineNum, kmer);
      lineNum += 1;
      continue;
    }
    count = strtol(countStr, &err, 10);
    if(*err != 0) {
      fprintf(stderr,
              "\tWARNING in line %d\ncould not read kmer count. Skipping\n",
              lineNum);
      lineNum += 1;
      continue;
    }


    //Load into hash
    if(count >= opt.lower && count <= opt.upper) {
      *k = kh_put(kmer, h, kmer,&absent);
      kh_key(h, *k) = strdup(kmer);
      //kh_value(h,*k) = count;
      *kmerCount += 1;
    }
    else *rkCount += 1;

    lineNum += 1;

  }

  free(line);
}


void count_filterReads(khash_t(kmer) *h,
                       khint_t k,
                       opts opt,
                       unsigned long *numReads,
                       unsigned long *numFiltered) {
  kseq_t *seq;
  int l;
  unsigned long int totalSeq = 0;
  char *kmer = malloc((opt.kmerlen + 1)*sizeof(char));
  kmer[opt.kmerlen] = '\0';

  seq = kseq_init(opt.seqFP);
  while ((l = kseq_read(seq)) >= 0) {
    if (l == 0) continue;
    totalSeq += l;
    if (l <= 0) {
      fprintf(stderr,"\nError in sequence file.\n");
      *numReads =  -1;
      break;
    }
    if(count_queryRead(h, k, opt, seq->seq.s, l, kmer)) {
      printf("@%s\n%s\n+\n%s\n",seq->name.s,seq->seq.s,seq->qual.s);
    }
    //Add filtered counter
    else *numFiltered += 1;
    *numReads += 1;
  }
  kseq_destroy(seq);
  free(kmer);
}


int *count_queryRead(khash_t(kmer) *h,
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
    else {
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
