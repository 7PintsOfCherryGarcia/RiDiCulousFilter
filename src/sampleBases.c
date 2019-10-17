#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include "../khash.h"
#include "../kseq.h"
#include "../ksort.h"

/*Initialize Heng Li's kseq
  KSEQ_INIT - Initialize sequence reading stream
  gzFile - Input file type
  gzread - Input file read function
*/
KSEQ_INIT(gzFile, gzread)

/*Initialize Heng Li's khash
  KHASH_SET_INIT_INT - Initiate hash table with ints as keys
  *DB - Name of hash table
  */
KHASH_SET_INIT_INT(intDB)

// Struct to store a sequence random index and lenght
typedef struct {
  unsigned int seqIdx;
  unsigned int seqLen;
  unsigned int seqRan;
} sequenceVal;

/*Initialize Heng Li's Ksort
*/
// Comparison function (less than)
#define pair_lt(a, b) ((a).seqRan < (b).seqRan)

//#define pair_lt(a, b) ((a).key < (b).key)
// Initialize with type(int) and comparison function(pair_lt)
KSORT_INIT(pair, sequenceVal, pair_lt)

// Don't know what this does
KSORT_INIT_GENERIC(int)

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
  sequenceVal *rnums;
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
  khash_t(intDB) *inth;
  inth = kh_init(intDB);
  khint_t intk;
  intk = kh_end(inth);
  int absent = 0;
  absent += 1;

  seq = kseq_init(seqFP);
  TotalSeq = 0;
  numSeq = 0;
  int maxNumSeq = 100000;
  //block controls reallocation size in case more space is needed
  int block = 2;
  rnums = (sequenceVal *)malloc(maxNumSeq*sizeof(sequenceVal));
  if(!rnums) {
    fprintf(stderr,"FAILED ALLOCATION\n");
    //TODO exit properly
    exit(-1);
  }
  int tmp = 0;
  while ((l = kseq_read(seq)) >= 0) {
    if (l == 0) continue;
    if (l <= 0) {
      fprintf(stderr,"\nError in sequence file.\n");
      break;
    }
    TotalSeq += l;
    //TODO bug if two sequences are assigend the same random number
    //TODO realloc if buffer runs out
    randInt = rand();
    if(randInt > tmp) {
      tmp = randInt;
    }
    //Reallocate more space if needed
    if(numSeq >= maxNumSeq && (numSeq%maxNumSeq) == 0) {
      fprintf(stderr,"Reallocating\n");
      rnums = realloc(rnums,(block*maxNumSeq)*sizeof(sequenceVal));
      block += 1;
    }

    rnums[numSeq].seqRan = randInt;
    rnums[numSeq].seqIdx = numSeq;
    rnums[numSeq].seqLen = l;
    intk = kh_put(intDB, inth, numSeq, &absent);
    numSeq += 1;
  }
  rnums = realloc(rnums,numSeq*sizeof(sequenceVal));
  fprintf(stderr,"Finished first pass\n");

  //Sort key array
  ks_mergesort(pair, numSeq, rnums, 0);

  //Loop over rnums and compute length until desired cummulative length
  unsigned int sampledLength = 0;
  unsigned int seqIdx = 0;
  while(sampledLength < 6000000) {
    sampledLength += rnums[seqIdx].seqLen;
      seqIdx += 1;
  }
  fprintf(stderr,"sampled length will be: %d\n",sampledLength);


  fprintf(stderr,"hash size before deleting is: %d\n",kh_size(inth));
  unsigned int seqKey;
  for(unsigned int i = seqIdx; i < numSeq; i++) {
    seqKey = rnums[i].seqIdx;
    intk = kh_get(intDB, inth, seqKey);
    kh_del(intDB,inth,intk);
  }
  fprintf(stderr,"hash size after deletingis: %d\n",kh_size(inth));

  numSeq = 0;
  seqFP = gzopen(seqfile,"r");
  seq = kseq_init(seqFP);
  unsigned int tmpLength = 0;
  while ((l = kseq_read(seq)) >= 0) {
    if (l == 0) continue;
    if (l <= 0) {
      fprintf(stderr,"\nError in sequence file.\n");
      break;
    }
    intk = kh_get(intDB,inth,numSeq);
    if (intk == kh_end(inth)) {
      ;;
    }
    else {
      tmpLength += l;
    }
    numSeq += 1;
  }
  fprintf(stderr,"Sampled number of seqs:%ld\n",numSeq);
  fprintf(stderr,"Sampled length: %d\n",tmpLength);
  fprintf(stderr,"Finished second pass\n");

  fprintf(stderr,"idx is: %d\n",seqIdx);
  fprintf(stderr, "%ld total bases in %ld sequences.\n", TotalSeq, numSeq);

  free(rnums);
  kseq_destroy(seq);
  kh_destroy(intDB, inth);
  gzclose(seqFP);
}
