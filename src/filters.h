

typedef struct {
  char kmer[99];
  int count;
} kmerCount;

typedef struct {
  int lower;
  int upper;
  char *file;
} opts;

int mainCountFilter(int argc, char *argv[]);

kmerCount *readKmerCounts(FILE *, kmerCount*);
