//zlib is used to read compressed data such as fastaFile.fq.gz
#include "zlib.h"
#include "khash.h"
KSEQ_INIT(gzFile, gzread)
/*Options structure
  Stores options to use for count command
*/

typedef struct {
  int minLen;           // lower kmer count bound (inclusive)
  int maxLen;           // upper kmer count bound (inclusive)
  char *seqfile;       // sequence filename
  gzFile seqFP;        // sequence filen object
} LENopts;

int main_length(int argc, char **argv);

void length_readOpt(int argc, char **argv, LENopts* opt);

void length_printOpt(LENopts opt);
