# HiChIP Data QC 
## Description
This is the description of the scripts that will perform QC steps in HiChIP data.

## Requirements

This script depends on the following tools in addition to the tools required for the alignment QC:

- [pysam](https://pysam.readthedocs.io/en/latest/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
- [deeptools](https://deeptools.readthedocs.io/en/develop/)
- [matplotlib](https://matplotlib.org/)
- [pandas](https://pandas.pydata.org/pandas-docs/stable/dsintro.html)

If you have already created conda environment using `create.sh` script in the `conda` folder, you have all these dependencies!

## Running


```
Usage: 
./hichip_qc.bash <reference_fasta> <read1_fastq> <reaf2_fastq>  <chipseq_peaks>  <output_prefix> <num_cores>

Example:
./hichip_qc.bash reference.fasta read1.fastq.gz reads2.fastq.gz chipseq_peaks.bed NA12878 8
```

`chipseq_peaks.bed` is a list of peaks called using ChipSeq data. We use this data from Encode data portal. 

## Output
This will print output as follows: 

```
Mapping Quality Threshold         : 40
Read1                             : 70,958,590
Read2                             : 70,868,593
Mapped pairs                      : 70,958,590
PCR dupe pairs                    : 8,852,388
Mapped nondupe pairs              : 62,106,202
Valid Pairs (cis>1000bp + trans)  : 30,988,708
Mapped nondupe pairs cis          : 51,114,356
Mapped nondupe pairs cis <=1000bp : 31,117,494
Mapped nondupe pairs cis >1000bp  : 19,996,862
Mapped nondupe pairs cis >10000bp : 6,399,456
Mapped nondupe trans pairs        : 10,991,846
Expected unique pairs at 300M sequencing:  217914791.5
Total ChIP peaks:       38,947
Mean ChIP peak size:    305 bp
Median ChIP peak size:  350 bp
Total reads in blacklist regions:       55,612(0.04%)
Total reads  in peaks:  7,746,420(5.33%)
Total reads in 500 bp around peaks:     14,328,101(9.85%)
Total reads in 1000 bp around peaks:    19,221,132(13.22%)
Total reads in 2000 bp around peaks:    27,841,150(19.15%)
Observed/Expected ratio for reads in peaks:     13.82
Observed/Expected ratio for reads in 500bp around  peaks:       9.7
Observed/Expected ratio for reads in 1000bp around  peaks:      8.03
Observed/Expected ratio for reads in 2000bp around  peaks:      6.59
```

Along with these statistics, the QC pipeline will output two plots. The first one is for the coverage enrichent around ChIP peaks. It would look as follows


![ChIP Enrichment Plot ](chip_enrichment_plot.png)

The pipeline also generates bigwig file for coverage tracks inferred from the BAM file in the output directory. 
