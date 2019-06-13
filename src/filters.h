#include "zlib.h"
#include "khash.h"

//Initialize Heng Li's khash and kseq
KHASH_MAP_INIT_STR(kmer, unsigned int)
KSEQ_INIT(gzFile, gzread)

typedef struct {
  char kmer[99];
  int count;
} kmerCount;


typedef struct {
  int kmerlen;
  int lower;
  int upper;
  char *kmerfile;
  char *seqfile;
  gzFile seqFP;
  double minfraction;
} opts;


int main_CountFilter(int, char **);


int readKmerCounts(FILE*, int, khash_t(kmer)*, khint_t*);


unsigned int loadKmerHash(khash_t(kmer)*, khint_t, kmerCount*, int);


unsigned long filterReads(khash_t(kmer)*, khint_t, opts);

int *query_read(khash_t(kmer)*,
                khint_t,
                opts,
                char*,
                int,
                char*);


void hash_print(khash_t(kmer)*, khint_t);


void hash_destroy(khash_t(kmer) *, khint_t);
