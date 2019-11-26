# Omni-C QC
## Description

This repo hosts a shell script (`omni-c_qc.bash`), that can be used to perform quick QC on shallowly sequenced Omni-C libraries.

## Requirements

This script depends on the following tools:

- [BWA](https://github.com/lh3/bwa)
- [samtools](https://github.com/samtools)
- [samblaster](https://github.com/GregoryFaust/samblaster)
- [preseq](http://smithlabresearch.org/software/preseq/)

All of these tools are installable via [BioConda](https://bioconda.github.io).

Install them however is most convinent for you. They are expected to be in your path.

## Running
Given paired FASTQ's and a reference FASTA, run:

```
./omni-c_qc.bash reference.fasta read1.fastq.gz read2.fastq.gz output.bam READGROUP_NAME
```

Substitute appropriate file names. READGROUP name is arbitrary.

## Output

After the script completes, it will print:

```
Read1                             : 5417421
Read2                             : 5417421
Mapped pairs                      : 5400976
PCR dupe pairs                    : 12293 0.00
Mapped nondupe pairs              : 5388683
Valid Pairs (cis>1000bp + trans)  : 4930332 0.91
Mapped nondupe pairs cis          : 4090918 0.76
Mapped nondupe pairs cis <=1000bp : 458351 0.09
Mapped nondupe pairs cis >1000bp  : 3632567 0.67
Mapped nondupe pairs cis >10000bp : 2751079 0.51
Mapped nondupe trans pairs        : 1297765 0.24
Expected unique pairs at 300M sequencing:  268938669.0
```

We consider a library to be acceptable if:

- cis pairs >1000 bp is greater than 20% of the total cis pairs.
- Expected unique pairs at 300M sequencing is ~120 million
