# Omni-C QC

## Description

This repo hosts a shell script (`qc.bash`), that can be used to perform quick QC on shallowly sequenced Omni-C libraries.

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
./qc.bash reference.fasta read1.fastq.gz read2.fastq.gz output.bam READGROUP_NAME
```

Substitute appropriate file names. READGROUP name is arbitrary.

## Output

After the script completes, it will print:

```
valid pairs          : 1728579
cis pairs            : 1891025
cis pairs >1000 bp   : 1092159
cis pairs >10000 bp  : 978995
trans pairs          : 636420

Expected unique pairs at 300M sequencing:  229193577.2
```

We consider a library to be acceptable if:

- cis pairs >1000 bp is greater than 20% of the total (usually much higher)
- Expected unique pairs at 300M sequencing is ~120 million
