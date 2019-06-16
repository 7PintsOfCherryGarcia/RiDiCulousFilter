#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include "kseq.h"
#include "main_gc.h"
#include "common.h"


/*
Usage message in case user specifies -h flag or ran with wrong options.
For example: -g 1.25. GC content can not exceed 100%
*/
int gc_usage() {
  fprintf(stderr,
          "Usage: RiDiCulous gc [options] -f SEQ_FILE\n");
  fprintf(stderr,"Options:\n\n");
  fprintf(stderr,"\t-g\tminimum GC content. Sequences must have at least \n");
  fprintf(stderr,"\t  \tthis GC content.\n\n");
  fprintf(stderr,"\t-f\tSequence file. fasta/q file with sequences or\n");
  fprintf(stderr,"\t  \tsequencing reads. Can be gziped. provide \"-\" if\n");
  fprintf(stderr,"\t  \treading from stdin. If using sequencing data eg.\n");
  fprintf(stderr,"\t  \tillumina reads, provide the -C flag to use canonical\n");
  fprintf(stderr,"\t  \tform.\n\n");
  fprintf(stderr,"\t-G\tMaximum GC content. Sequences may have at most.\n");
  fprintf(stderr,"\t  \tthis GC content.\n\n");
  fprintf(stderr,"\t-C\tUse canonical kmers. Sequences are reversed complemente");
  fprintf(stderr,"d\n\t \tand smallest lexicographically sequence is used\n\n");
  fprintf(stderr,"\t-h\tThis help message\n");
  fprintf(stderr,"\n");
  return -1;
}


void gc_readOpt(int argc, char **argv, GCopts *opt) {
  int elem;
  while (( elem = getopt(argc, argv, "Cg:G:f:h") ) >= 0) {
    switch(elem) {
    case 'f':
      opt->seqfile = optarg;
      break;
    case 'h':
      exit(gc_usage());
    case 'g':
      opt->minGC = atof(optarg);
      if(!opt->minGC || opt->minGC < 0 || opt -> minGC > 1) {
        fprintf(stderr,"\t ERROR: Option \"-g %s\" not valid\n\n",optarg);
        exit(gc_usage());
      }
      break;
    case 'G':
      opt->maxGC = atof(optarg);
      if(!opt->maxGC || opt->maxGC < 0 || opt -> maxGC > 1) {
        fprintf(stderr,"\t ERROR: Option \"-G %s\" not valid\n\n",optarg);
        exit(gc_usage());
      }
      break;
    case 'C':
      opt->canonical=true;
    }
  }

  if (!opt->seqfile) {
    fprintf(stderr,"\t ERROR: Please a sequence file (-f)\n\n");
    exit(gc_usage());
  }

  if(opt->minGC > opt->maxGC) {
    fprintf(stderr, "\tERROR: maximum GC content must be greater than");
    fprintf(stderr," minimum gc content.\n\n");
    exit(gc_usage());
  }


  gc_printOpt(*opt);
}


/*
  Print input options
*/
void gc_printOpt(GCopts opt) {
  fprintf(stderr,"\t sequence file: %s\n",opt.seqfile);
  fprintf(stderr,"\t minimum GC content: %f\n",opt.minGC);
  fprintf(stderr,"\t maximum GC content: %f\n",opt.maxGC);
  if(opt.canonical) {
    fprintf(stderr,"\t use canonical kmers: true\n");
  }
  else {
    fprintf(stderr,"\t use canonical kmers: false\n");
  }
}


int main_gc(int argc, char **argv) {
  GCopts opt;
  opt.seqfile = NULL;
  opt.seqFP = NULL;
  opt.minGC = 0.40;
  opt.maxGC = 0.60;
  opt.canonical = false;

  //Read options
  gc_readOpt(argc, argv, &opt);

  //Open sequence file
  if(strcmp(opt.seqfile,"-") == 0) {
    opt.seqfile = "/dev/stdin";
  }

  opt.seqFP = gzopen(opt.seqfile,"r");
  if(!opt.seqFP) {
    fprintf(stderr,"\t ERROR: Can't open sequence file.\n\n");
    return gc_usage();
  }

  //Loop over reads apply filter and write to stdout
  fprintf(stderr,"\t Filtering reads.\n");
  unsigned long numReads = 0;
  unsigned long numFiltered = 0;
  gc_filterReads(opt, &numReads, &numFiltered);
  fprintf(stderr,"\t%ld reads processed\n",numReads);
  fprintf(stderr,"\t%ld reads filtered\n",numFiltered);
  fprintf(stderr,"\t%ld reads kept\n",numReads-numFiltered);
  gzclose(opt.seqFP);
  return 0;
}


void gc_filterReads(GCopts opt,
                    unsigned long *numReads,
                    unsigned long *numFiltered) {
  kseq_t *seq;
  int l;
  unsigned long int totalSeq = 0;
  char *revComp = NULL;

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

    if(gc_queryRead(opt, seq->seq.s, l)) {
      //TODO Correct for qhen input is fasta
      printf("@%s\n%s\n+\n%s\n",seq->name.s,seq->seq.s,seq->qual.s);
    }
    //Add filtered counter
    else *numFiltered += 1;
    *numReads += 1;
  }
  free(revComp);
  kseq_destroy(seq);
}


int *gc_queryRead(GCopts opt, char *seq,int l) {
  double GCFraction = 0;
  int GCcounter = 0;
  //Loop over bases
  for(int i = 0; i < l; i++) {
    if(*&seq[i] == 71 || *&seq[i] == 67 || *&seq[i] == 103 || *&seq[i] == 99) {
      GCcounter += 1;
    }
  }
  GCFraction = (double) GCcounter / (double) l;
  if (GCFraction >= opt.minGC && GCFraction <= opt.maxGC) {
    return &true;
  }
  else {
    return NULL;
  }
}
