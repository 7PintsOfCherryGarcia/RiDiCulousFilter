

typedef struct {
  char *kmer;
  int count;
} kmerCount;

int readKmerCounts(FILE *, kmerCount*);
