#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "khash.h"
#include "kseq.h"
#include "main_length.h"
#include "common.h"

int length_usage() {
  fprintf(stderr,
          "Usage: RiDiCulous length [options] -f SEQ_FILE\n");
  fprintf(stderr,"Options:\n\n");
  fprintf(stderr,"\t-f\tSequence file. fasta/q file with sequences or\n");
  fprintf(stderr,"\t  \tsequencing reads. Can be gziped. provide \"-\" if\n");
  fprintf(stderr,"\t  \treading from stdin. If using sequencing data eg.\n");
  fprintf(stderr,"\t  \tillumina reads, provide the -C flag to use canonical\n");
  fprintf(stderr,"\t  \tform.\n\n");
  fprintf(stderr,"\t-l\tMinimum sequence length. Only sequences longer or \n");
  fprintf(stderr,"\t  \tequal to this value will be kept.\n\n");
  fprintf(stderr,"\t-u\tMaximum sequence length. Only sequences shorter or\n");
  fprintf(stderr,"\t  \tequal to thie value will be kept\n\n");
  fprintf(stderr,"\t-h\tThis help message\n");
  fprintf(stderr,"\n");
  return -1;
}


void length_readOpt(int argc, char **argv, LENopts* opt) {
  int elem;
  while (( elem = getopt(argc, argv, "f:hk:l:m:u:") ) >= 0) {
    switch(elem) {
    case 'l':
      opt->minLen = atoi(optarg);
      if(!opt->minLen) {
        fprintf(stderr,"\t ERROR: Option \"-l %s\" not valid\n\n",optarg);
        exit(length_usage());
      }
      break;
    case 'u':
      opt->maxLen = atoi(optarg);
      if(!opt->maxLen) {
        fprintf(stderr,"\t ERROR: Option \"-u %s\" not valid\n\n",optarg);
        exit(length_usage());
      }
      break;
    case 'f':
      opt->seqfile = optarg;
      break;
    case 'h':
      exit(length_usage());
      break;
    }
  }

  if (!opt->seqfile) {
    fprintf(stderr,"\t ERROR: Please provide a sequence file (-f)\n\n");
    exit(length_usage());
  }

  if((opt->minLen < 0)) {
    fprintf(stderr, "\tERROR: please provide minimum sequence length (-l) ");
    fprintf(stderr, "and/or maximum sequence length (-u).\n\n");
    exit(length_usage());
  }

  if(!opt->maxLen) {
    opt->maxLen = 10000000000;
  }
  else if( opt->minLen >= opt->maxLen) {
    fprintf(stderr, "\tERROR: minumum seq length must be less than maximum length.\n\n");
    exit(length_usage());
  }



  length_printOpt(*opt);
}


/*
  Print input options
*/
void length_printOpt(LENopts opt) {
  fprintf(stderr,"\t sequence file: %s\n",opt.seqfile);
  fprintf(stderr,"\t minimum sequence length: %ld\n",opt.minLen);
  fprintf(stderr,"\t maximum sequence length: %ld\n",opt.maxLen);
}

int main_length(int argc, char **argv) {
  //Set options
  LENopts opt;
  opt.minLen = -1;
  opt.maxLen  = 0;
  opt.seqfile = NULL;
  opt.seqFP = NULL;

  //Read options
  length_readOpt(argc - 1, argv + 1, &opt);

  //Open sequence file
  if(strcmp(opt.seqfile,"-") == 0) {
    opt.seqfile = "/dev/stdin";
  }

  opt.seqFP = gzopen(opt.seqfile,"r");
  if(!opt.seqFP) {
    fprintf(stderr,"\t ERROR: Can't open sequence file.\n\n");
    return length_usage();
  }

  kseq_t *seq;
  int l;
  unsigned long int totalSeq = 0;
  unsigned long int numSeq = 0;
  seq = kseq_init(opt.seqFP);
  while ((l = kseq_read(seq)) >= 0) {
    if (l == 0) continue;

    if(l >= opt.minLen && l <= opt.maxLen) {
      printf(">%s\n%s\n",seq->name.s,seq->seq.s);
    }

    totalSeq += l;
    if (l <= 0) {
      fprintf(stderr,"\nError in sequence file.\n");
      break;
    }
    numSeq += 1;

  }
  return 0;
}
