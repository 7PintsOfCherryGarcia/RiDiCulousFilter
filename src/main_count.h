//zlib is used to read compressed data such as fastaFile.fq.gz
#include "zlib.h"
//Heng Li's hash table implementation
#include "khash.h"


/*Initialize Heng Li's khash
  KHASH_MAP_INIT_STR - Initiate hash table with string pointers as keys
  kmer - Name of hash table
  unsigend int - type of values to use
 */
KHASH_MAP_INIT_STR(kmer, unsigned int)

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
  int kmerlen;         // kmer length
  int lower;           // lower kmer count bound (inclusive)
  int upper;           // upper kmer count bound (inclusive)
  char *kmerfile;      // kmer count filename
  char *seqfile;       // sequence filename
  gzFile seqFP;        // sequence filen object
  double minfraction;  // kmer coverage fraction threshold
  int canonical;       // If set to 1, canonical kmers are extracted
  int kmerCMD;         // If set to 1, kmer comand is used instead of count
} COUNTopts;

/*
Main function for count and kmer commands
INPUT:
    argc - int, number of command line arguments
    argv - char**, string of each command line argument
RETURNS:
    int, 0 on success, -1 on failure
*/
int main_count(int argc, char **argv);


/*
Parses command arguments, sets appropiate values and checks for correct input
INPUT:
    argc - int, number of command line arguments
    argv - char**, string of each command line argument
    opt  - COUNTopts*, pointer to COUNTopts struct containing parameters
RETURNS:
    nothing
*/
void count_readOpt(int, char**, COUNTopts*);


/*
*/
void count_printOpt(COUNTopts);


/*
*/
void count_readKmers(FILE*, COUNTopts, khash_t(kmer)*, khint_t*, int*, int*);


/*
*/
void count_filterReads(khash_t(kmer)*,
                       khint_t,
                       COUNTopts,
                       unsigned long*,
                       unsigned long*);


/*
*/
int *count_queryRead(khash_t(kmer)*,
                     khint_t,
                     COUNTopts,
                     char*,
                     int,
                     char*);


/*
*/
void hash_print(khash_t(kmer)*, khint_t);


/*
*/
void hash_destroy(khash_t(kmer) *, khint_t);
