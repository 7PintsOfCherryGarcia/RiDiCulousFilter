# RiDiCulousFilter
## Filter sequences based on a variety of conditions

## About
RiDiCulousFilter lets you select nucleotide sequence data such as sequencing reads or assembly contigs
based on different criteria. For example: kmer count content or GC content. In the future 
other filtering parameters will be added such as protein domain presence, and sequence 
signatures. (I want reads/sequences that start in AT and end in G for example). For now
only GC content and kmer count/presense are implemented.
The backbone of RiDiCulousFilter is the wonderfull [klib](https://github.com/attractivechaos/klib) library by [Attractive Chaos](https://github.com/attractivechaos)

## How does it work
It's quite simple. RiDiCulousFilter loops over all your sequences/reads, 
computes the desired filter and determines if conditions are met. If so,
sequences are printed to stdout. You'll have a better understandig if you
dig into the code I have tried to document it as much as possible.

## How do I get this thing
Easy peasy lemon squeezy. You will need any linux machine (should work with IOS).
Other requirements are the zlib (to read compressed input) library, gcc (to compile), and cmake (generate the makefile).
RiDiCulous also uses khash and kseq from [klib](https://github.com/attractivechaos/klib). These .h files are already included in this repo. 

```
git clone https://github.com/7PintsOfCherryGarcia/RiDiCulousFilter.git
cd RiDiCulousFilter
cmake .
make
```

## How do I use RiDiCulousFilter?

You will find it is straight forward and simple. Here are a couple of examples.

### Filter sequencing reads based on kmer counts
Imagine a fastq file: [Microverse](https://rickandmorty.fandom.com/wiki/Microverse).fq.gz. It contains reads of a whole genome shotgun sequencing experiment of an indivdual of the species [*Zeep Zanflorp*](https://rickandmorty.fandom.com/wiki/Zeep_Xanflorp), an interesting organism with a 10Gb genome. We know Microverse.fq.gz contains data such that *Zeep Zanflorp* was sequenced at 30X depth. We are interested in genes at high copy number say 8, known to be involved in the synthesis of [Gooble Boxes](https://rickandmorty.fandom.com/wiki/Gooble_Box). Instead of trying to assemble the entire genome and then search for such high copy genes, we extract all reads that have kmers at a count expected from regions at 8 copies.

**Determining kmer count threshold***

8(copies)*30X = 240, ***Reads from these regions should have kmers whose counts are around 240***

***Count kmers. We will use [kmc](https://github.com/refresh-bio/KMC) and a kmer length of 25***

```
kmc -k25 -cs500 Microverse.fq.gz Microverse
```

***Run RiDiCulousFilter on the output kmer count table of kmc***

```
kmc_dump Microverse /dev/stdout | RiDiCulous count -c - -f Microverse.fq.gz -l 220 -u 260 -C -m 0.8 -k 25 > filteredMicroverse.fq
```

filteredMicroverse.fq will contain all the reads that had at least 80% of it's kmers at a count between 220 and 260 inclusive


# FAQ and comments

> Hey, I don't want/can't to use cmake!!!

Fine, just compile with

```
gcc -Wall -O3 -o RiDiCulous src/*.c -lm -I. -lz 
```

> Bro, I need to filter my sequences based on "some outreagous filter that nobody will never ever use". 

Sure Fam!! No problemo just hit me up in the [issues](https://github.com/7PintsOfCherryGarcia/RiDiCulousFilter/issues) section. 

> Your "How does it work" section sucks. Please elaborate!!

Sorry to hear that homie. I have put significant effort in documenting almost every single line of code.
Try to dig in, you will it easy, I promise :)
