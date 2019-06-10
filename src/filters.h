#include "khash.h"

//Initialize Heng Li's khash
KHASH_MAP_INIT_STR(kmer, unsigned int)

typedef struct {
  char kmer[99];
  int count;
} kmerCount;

typedef struct {
  int lower;
  int upper;
  char *file;
} opts;

int mainCountFilter(int, char **);

kmerCount *readKmerCounts(FILE*, kmerCount*, int*);

int loadKmerHash(khash_t(kmer)*, khint_t, kmerCount*, int);
