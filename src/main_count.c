#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "khash.h"
#include "kseq.h"
#include "main_count.h"
#include "common.h"



/*
Usage message in case user specifies -h flag or ran with wrong options.
For example: -l 10 -u 4. Upper bound must be greater than lower bound
*/
int count_usage() {
  fprintf(stderr,
          "Usage: RiDiCulous count [options] -c KMER_FILE -f SEQ_FILE\n");
  fprintf(stderr,"Options:\n\n");
  fprintf(stderr,"\t-c\tkmer count file. Tab separated text file with\n");
  fprintf(stderr,"\t  \tkmer sequences as first column and counts in\n");
  fprintf(stderr,"\t  \tsecond column.\n\n");
  fprintf(stderr,"\t-C\tUse canonical kmers. Sequences are reversed compliment");
  fprintf(stderr,"ed\n\t \tand smallest lexicographically sequence is used\n\n");
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
  fprintf(stderr,"\t  \tless than or equal to this value will be consife");
  fprintf(stderr,"red.\n\n");
  fprintf(stderr,"\t-h\tThis help message\n");
  fprintf(stderr,"\n");
  return -1;
}


/*
Parse options for "count" command. See main_count.h file.
*/
void count_readOpt(int argc, char **argv, COUNTopts* opt) {
  int elem;
  while (( elem = getopt(argc, argv, ":c:Cf:hk:l:m:u:") ) >= 0) {
    switch(elem) {
    case 'l':
      opt->lower = atoi(optarg);
      if(!opt->lower) {
        fprintf(stderr,"\t ERROR: Option \"-l %s\" not valid\n\n",optarg);
        exit(count_usage());
      }
      break;
    case 'u':
      opt->upper = atoi(optarg);
      if(!opt->upper) {
        fprintf(stderr,"\t ERROR: Option \"-u %s\" not valid\n\n",optarg);
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
        fprintf(stderr,"\t ERROR: Option \"-m %s\" not valid\n\n",optarg);
        exit(count_usage());
      }
      break;
    case 'k':
      opt->kmerlen = atoi(optarg);
      if(!opt->kmerlen) {
        fprintf(stderr,"\t ERROR: Option \"-k %s\" not valid\n\n",optarg);
        exit(count_usage());
      }
      break;
    case 'C':
      opt->canonical=true;
    }
  }

  if (!opt->kmerfile || !opt->seqfile) {
    fprintf(stderr,"\t ERROR: Please provide kmer count file (-c) and\n");
    fprintf(stderr,"\t        sequence file (-f)\n\n");
    exit(count_usage());
  }

  if((opt->lower < 0) || (opt->upper < 0 )) {
    fprintf(stderr, "\tERROR: please provide lower kmer count bound (-l) ");
    fprintf(stderr, "and upper kmer count bound (-u).\n\n");
    exit(count_usage());
  }

  if( opt->lower >= opt->upper) {
    fprintf(stderr, "\tERROR: lower bound must be less than upper bound.\n\n");
    exit(count_usage());
  }


  if(opt->minfraction <= 0 || opt->minfraction > 1) {
    fprintf(stderr,"\t WARNING: minimum fraction not in \(0,1] range.\n");
    fprintf(stderr,"\t         Setting to default -m 0.8\n\n");
    opt->minfraction = 0.80;
  }

  count_printOpt(*opt);
}


/*
Print input options
*/
void count_printOpt(COUNTopts opt) {
  fprintf(stderr,"\t kmer count file: %s\n",opt.kmerfile);
  fprintf(stderr,"\t sequence file: %s\n",opt.seqfile);
  fprintf(stderr,"\t lower count bound: %d\n",opt.lower);
  fprintf(stderr,"\t upper count bound: %d\n",opt.upper);
  fprintf(stderr,"\t minimum read coverage fraction: %f\n",opt.minfraction);
  fprintf(stderr,"\t kmer length: %d\n",opt.kmerlen);
  if(opt.canonical) {
    fprintf(stderr,"\t use canonical kmers: true\n");
  }
  else {
    fprintf(stderr,"\t use canonical kmers: false\n\n");
  }
}


/*
Main function for command "count". See main_count.h file.
*/
int main_count(int argc, char **argv) {
  //Set options
  COUNTopts opt;
  opt.kmerfile = NULL;
  opt.seqfile = NULL;
  opt.seqFP = NULL;
  opt.lower = -1;
  opt.upper = -1;
  opt.minfraction = 0.80;
  opt.kmerlen = 31;
  opt.canonical = false;

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
    fprintf(stderr,"\t ERROR: Can't open input file %s.\n\n", opt.kmerfile);
    return(count_usage());
  }

  //Open sequence file
  if(strcmp(opt.seqfile,"-") == 0) {
    if(stdinFlag == 1) {
      fprintf(stderr,"\t ERROR: -f - and  -c - are incompatibe.\n");
      fprintf(stderr,"\t Only one input file can be read from stdin.\n\n");
      return count_usage();
    }
    opt.seqfile = "/dev/stdin";
  }

  opt.seqFP = gzopen(opt.seqfile,"r");
  if(!opt.seqFP) {
    fprintf(stderr,"\t ERROR: Can't open sequence file.\n\n");
    return count_usage();
  }



  //Declare hash objects
  khash_t(kmer) *h;
  h = kh_init(kmer);
  khint_t k;
  k = kh_end(h);

  //Load kmer counts into hash table
  fprintf(stderr,"\t Loading kmer count table\n");
  int kmerCount = 0;
  int rkCount = 0;
  count_readKmers(fp, opt, h, &k, &kmerCount, &rkCount);
  if(kmerCount <= 0) {
    fprintf(stderr, "\t ERROR: kmer count table could not be parsed.\n\n");
    return count_usage();
  }
  fclose(fp);
  fprintf(stderr,"\t\t Loaded %d kmers.\n",kmerCount);
  fprintf(stderr,"\t\t Skipped %d kmers.\n",rkCount);


  //Loop over reads apply filter and write to stdout
  fprintf(stderr,"\t Filtering reads.\n");
  unsigned long numReads = 0;
  unsigned long numFiltered = 0;
  count_filterReads(h, k, opt, &numReads, &numFiltered);
  fprintf(stderr,"\t%ld reads processed\n",numReads);
  fprintf(stderr,"\t%ld reads filtered\n",numFiltered);
  fprintf(stderr,"\t%ld reads kept\n",numReads-numFiltered);

  //Free memoory
  hash_destroy(h, k);
  gzclose(opt.seqFP);
  return 0;
}


/*
Parse kmer and their counts. Store in hash table kmers where counts satisfy
lower and uper bound thersholds. See main_count.h file.
*/
void count_readKmers(FILE *fp,
                     COUNTopts opt,
                     khash_t(kmer) *h,
                     khint_t *k,
                     int *kmerCount,
                     int *rkCount) {

  //Pointer to line string read from kmer count file
  char *line = NULL;
  size_t n = 0;

  //Pointer to kmer string extracted from first part of line
  char *kmer;
  //Pointer to count string extracted from second part of line
  char *countStr;
  //Count variable (value represented in *countStr is stored here)
  int count;
  //Error variable when for when parsing goes wrong
  char *err;
  //Number of characters read in each line. Used for removing newline character
  int numChars;
  //Keep track of line number for error reporting
  int lineNum = 1;
  int absent = 0;

  while((numChars = getline(&line, &n, fp)) != -1) {
    line[numChars - 1] = '\0';
    kmer = strtok(line, "\t");
    if(!kmer) {
      fprintf(stderr,"\t\t WARNING: In line %d\n\t\t\t Could not read", lineNum);
      fprintf(stderr," kmer sequence. Skipping\n");
      lineNum += 1;
      continue;
    }
    else if(strlen(kmer) != opt.kmerlen) {
      fprintf(stderr,"\t\t WARNING in line %d\n\t\t\t kmer length of",lineNum);
      fprintf(stderr," different size. Skipping\n");
      lineNum += 1;
      continue;
    }
    countStr = strtok(NULL,"\t");
    if(!countStr) {
      fprintf(stderr,"\t\t WARNING: In line %d\n\t\t\t kmer %s", lineNum, kmer);
      fprintf(stderr," without a count. Skipping\n");
      lineNum += 1;
      continue;
    }
    count = strtol(countStr, &err, 10);
    if(*err != 0) {
      fprintf(stderr,"\t\t WARNING in line %d\n\t\t\tcould not", lineNum);
      fprintf(stderr," read kmer count. Skipping\n");
      lineNum += 1;
      continue;
    }


    //Load into hash
    if(count >= opt.lower && count <= opt.upper) {
      *k = kh_put(kmer, h, kmer,&absent);
      kh_key(h, *k) = strdup(kmer);
      *kmerCount += 1;
    }
    else *rkCount += 1;

    lineNum += 1;

  }

  free(line);
}


void count_filterReads(khash_t(kmer) *h,
                       khint_t k,
                       COUNTopts opt,
                       unsigned long *numReads,
                       unsigned long *numFiltered) {
  kseq_t *seq;
  int l;
  unsigned long int totalSeq = 0;
  char *kmer = malloc((opt.kmerlen + 1)*sizeof(char));
  kmer[opt.kmerlen] = '\0';
  char *revComp;
  seq = kseq_init(opt.seqFP);
  while ((l = kseq_read(seq)) >= 0) {
    if (l == 0) continue;
    totalSeq += l;
    if (l <= 0) {
      fprintf(stderr,"\nError in sequence file.\n");
      *numReads =  -1;
      break;
    }

    if(opt.canonical) {
      revComp = gc_revComp(l, seq->seq.s);
      //Compare sequence and reverse compliment. Keep lexicgraphically
      //smales one
      seq->seq.s = (compSeq(seq->seq.s,revComp,l) <= 0) ? seq->seq.s:revComp;
    }

    if(count_queryRead(h, k, opt, seq->seq.s, l, kmer)) {
      printf("@%s\n%s\n+\n%s\n",seq->name.s,seq->seq.s,seq->qual.s);
    }
    //Add filtered counter
    else *numFiltered += 1;
    *numReads += 1;
  }
  free(revComp);
  kseq_destroy(seq);
  free(kmer);
}


int *count_queryRead(khash_t(kmer) *h,
                     khint_t k,
                     COUNTopts opt,
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
    return &false;
  }
  return NULL;
}


//For testing purposes
void hash_print(khash_t(kmer) *h, khint_t k) {
  //Loop over hash integers and print key value pair
  for (k = 0; k < kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      printf("%s\t%d\n",kh_key(h,k),kh_val(h,k));
    }
  }
}


/*
Free string pointers used as keys in hash table. See main_count.h
*/
void hash_destroy(khash_t(kmer) *h, khint_t k) {
  for (k = 0; k < kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      free((char*)kh_key(h, k));
    }
  }
  kh_destroy(kmer, h);
}
