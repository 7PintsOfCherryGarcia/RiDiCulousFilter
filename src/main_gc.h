//zlib is used to read compressed data such as fastaFile.fq.gz
#include "zlib.h"

/*Initialize Heng Li's kseq
  KSEQ_INIT - Initialize sequence reading stream
  gzFile - Input file type
  gzread - Input file read function
*/
KSEQ_INIT(gzFile, gzread)

/*Options structure
  Stores options to use for count command
*/

typedef struct {
  double minGC;        // minimum GC content thershold
  double maxGC;        // maximum GC content thershold
  char *seqfile;       // sequence filename
  gzFile seqFP;        // sequence file object
} GCopts;



int main_gc(int, char **);


void gc_printOpt(GCopts);


void gc_readOpt(int, char**, GCopts*);


void gc_filterReads(GCopts, unsigned long*, unsigned long*);


int *gc_queryRead(GCopts, char*, int);



