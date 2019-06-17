# RiDiCulousFilter
##@ Filter sequences based on a variety of conditions (even if they are comical)

## About
RiDiCulousFilter lets you select sequence data such as sequencing reads or assembly contigs
based on different criteria. For example: kmer count content or GC content. In the future 
other filtering parameters will be added such as protein domain prescence, and sequence 
signatures. (I want reads/sequences that start in AT and end in G for example). For now
only GC content and kmer count/presense are implemented.
The backbone of RiDiCulousFilter is the wonderfull [klib](https://github.com/attractivechaos/klib) library by [Attractive Chaos](https://github.com/attractivechaos)

## How does it work
It's quite simple. RiDiCulousFilter loops over all your sequences/reads, 
computes the desired filter and determines if conditions are met. If so,
sequences are printed to stdout.

## How do I get this dumb thing
Easy peasy lemon squasy. You will need any linux machine (should work with IOS).
Other requirements are the zlib library, gcc, and cmake.

```
git clone https://github.com/7PintsOfCherryGarcia/RiDiCulousFilter.git
cd RiDiCulousFilter
cmake .
make
```

## How do I use RiDiCulousFilter?
You will find it is straight forward and simple.
### Filter sequencing reads based on kmer counts
Imagine we have a fastq file myLibrary.fq.gz which contains the sequencing reads of a whole genome shotgun sequencing experiment of an indivdual of the species [*Zeep Zanflorp*](https://rickandmorty.fandom.com/wiki/Zeep_Xanflorp), a very interesting organism with a 10Gb genome. We know myLibrary.fq.gz contains enough data such that *Zeep Zanflorp* was sequenced at 30X depth. For our purposes, we are interested genes at high copy number say 8, which are known to be involved in the synthesis of Gooble Boxes. Instead of trying to assemble the entire genome and then search for such high copy genes, we want to extract all sequencing reads 



# FAQ and comments

> Hey, I don't want to use your dumb cmake!!!

Fine, just compile with

```
gcc -Wall -O3 -o RiDiCulous src/*.c -lm -I. -lz 
```

> Bro, I need to filter my sequences based on "some outreagous filter that nobody will never ever use". 

Sure Fam!! No problemo just hit me up in the [issues](https://github.com/7PintsOfCherryGarcia/RiDiCulousFilter/issues) section. 

> Your "How does it work" section sucks. Please elaborate!!

Sorry to hear that homie. I have put significant effort in documenting almost every single line of code.
Try to dig in, you will it quite straight forward :)
