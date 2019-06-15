#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "khash.h"
#include "kseq.h"
#include "main_count.h"


//Used as a TRUE/FALSE flag
static int true = 1;


/*
Usage message in case user specifies -h flag or ran with wrong options.
For example: -c 1.25. GC content can not exceed 100%
*/
int count_usage() {
  fprintf(stderr,
          "Usage: RiDiCulous gc [options] -f SEQ_FILE\n");
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
