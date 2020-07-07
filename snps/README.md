
# Calling SNPs with OmniC data

Snps can be called with the included scripts like: 

```
snps_pipe.bash sample.bam outdir/rootname ConfidentRegions.bed.gz GRCh38.p12.fa 
```

Where 'rootname' is the root name for all the output files and ConfidentRegions.bed is a file describing the intervals over which SNPs should be called.  

This produces a raw set of snps and a separate high confidence set of snps.  High confident snps are snps that fall in regions where there are no MQ=0 aligned reads.  Some regions that are confident regions for shotgun are not confident regions for OmniC because the pair size distribution of *-C data is very broad and some mates whose placements could be disambiguated by bwa with shotgun data can not be disambiguated with *-C data.   A modified version of bwa, or some post-processing, might be able to constrain some of these mates based on the broader Omni-C pair size distribution but for now we just identify these difficult regions and exclude them from the results. 

  





