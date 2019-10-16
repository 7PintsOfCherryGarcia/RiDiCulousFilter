#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include "../khash.h"
#include "../kseq.h"

/*Initialize Heng Li's kseq
  KSEQ_INIT - Initialize sequence reading stream
  gzFile - Input file type
  gzread - Input file read function
*/
KSEQ_INIT(gzFile, gzread)

/*Initialize Heng Li's khash
  KHASH_MAP_INIT_STR - Initiate hash table with string pointers as keys
  kmer - Name of hash table
  unsigend int - type of values to use
*/
KHASH_MAP_INIT_INT(seqDB, unsigned int)
KHASH_MAP_INIT_INT(intDB, unsigned int)

int main(int argc, char **argv) {
  //Random seed
  //TODO Give option to not always being random
  srand(time(0)*time(0));


  char *seqfile;
  gzFile seqFP;
  unsigned long TotalSeq;
  unsigned long numSeq;
  kseq_t *seq;
  int l;
  int *rnums;
  unsigned int randInt;

  if(argc < 2) {
    fprintf(stderr,"ERROR\n");
  }

  for(int i = 0; i < argc; i++) {
    fprintf(stderr,"%s\n",argv[i]);
  }

  //Open seq file and initialize
  if(strcmp(argv[1],"-") == 0) {
    seqfile = "/dev/stdin";
  }
  else {
    seqfile = argv[1];
  }

  fprintf(stderr,"Seqfile: %s\n",seqfile);
  seqFP = gzopen(seqfile,"r");
  if(!seqFP) {
    fprintf(stderr,"ERROR\n");
    //TODO exit properly, free memory and destroy hash
    exit(-1);
  }

  //Start hash tables
  khash_t(seqDB) *h;
  khash_t(intDB) *inth;
  h = kh_init(seqDB);
  inth = kh_init(intDB);
  khint_t k;
  khint_t intk;
  k = kh_end(h);
  intk = kh_end(inth);
  int absent = 0;
  absent += 1;

  seq = kseq_init(seqFP);
  TotalSeq = 0;
  numSeq = 0;
  rnums = (int *)malloc(1000*sizeof(int));
  if(!rnums) {
    fprintf(stderr,"FAILED ALLOCATION\n");
    //TODO exit properly
    exit(-1);
  }
  while ((l = kseq_read(seq)) >= 0) {
    if (l == 0) continue;
    if (l <= 0) {
      fprintf(stderr,"\nError in sequence file.\n");
      break;
    }
    TotalSeq += l;
    numSeq += 1;
    //TODO bug if two sequences are assigend the same random number
    //TODO realloc if buffer runs out
    randInt = rand();
    rnums[numSeq] = randInt;
    intk = kh_put(intDB, inth, randInt, &absent);
    kh_key(inth, intk) = numSeq;
  }
  fprintf(stderr, "%ld total bases in %ld sequences.\n", TotalSeq, numSeq);
  fprintf(stderr,"%d\n",rnums[10]);




  //Loop over sequences, count number of sequences and total bases

  //k = kh_put(seqDB,h,s,&absent);




  //for (k = 0; k < kh_end(h); ++k) {
  //  if (kh_exist(h, k)) {
  //    free((char*)kh_key(h, k));
  //  }
  //}
  free(rnums);
  kseq_destroy(seq);
  kh_destroy(seqDB, h);
  kh_destroy(intDB, inth);
  gzclose(seqFP);
}
