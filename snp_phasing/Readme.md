# Omni-C SNV calling and phasing pipeline

In order to run this pipeline, you need to align and QC data using `omni-c_qc.bash` script from this repository. Once you have the BAM file, you can proceed with the SNP calling. 

## Environment

In order to run these scripts, you would require following dependencies:

- [tabix](https://anaconda.org/bioconda/tabix)
- [WhatsHap](https://whatshap.readthedocs.io/en/latest/)
- [HapCUT2](https://github.com/vibansal/HapCUT2)
- [java](https://anaconda.org/bioconda/java-jdk)
- [pysam](https://pysam.readthedocs.io/en/latest/index.html)

All the dependencies can be installed with Conda. Also, we assume you have all the depencies required for alignment/qc pipeline installed at this time. 

## SNV calling with Omni-C data

We use `HaplotypeCaller` program from GATK to call SNVs.

```
./snp_calling.sh <reference_fasta> <alignment_bam> <output_prefix>
```
This run will produce `vcf`, `vcf.gz`, and `vcf.gz.tbi` file in the end. 

If you want to filter this vcf based on `QD` and `FS` filters from the VCF INFO field, you can use `qd_fs_filter.py` script as follows:

```
./qd_fs_filter.py -qdthresh 2 -fsthresh 60 -vcfin in.vcf.gz -vcfout out.vcf
```
This would generate `out.vcf.gz` and `out.vcf.gz.tbi` files that are fitlered for given qd and fs thresholds. 

## Phasing variants with Omni-C data

To phase variants, we use HapCUT2 phasing program. The wrapper script to run for phasing is `phasing_workflow.sh`. It can be run as

```
./phasing_workflow.sh <variants_vcf> <alignment_bam> <output_prefix>
```


Note that `<vcf>` file here should NOT be gzipped as HapCUT2 requires unzipped VCF file. This script will first filter input VCF to contain just signel or biallelic variants as HapCUT2 operataes only on these types of variants. The output of the script is `<output_prefix>_hapcut_output.phased.VCF.gz` file that contains phased variants. 

## Evaluating SNPs and phasing for NA12878

If you want test data to run this pipeline, you can download the pre-aligned BAM (and its index) file from here: 

https://dovetail-public-data.s3-us-west-2.amazonaws.com/NA12878_data/NA12878_OmniC.bam


https://dovetail-public-data.s3-us-west-2.amazonaws.com/NA12878_data/NA12878_OmniC.bam.bai

After you download it, first run the SNP calling as follows as follows:
```
./snp_calling.sh hg38.fasta NA12878_OmniC.bam NA12878
```

This should produce `NA12878_variants.vcf`, `NA12878_variants.vcf.gz`, and `NA12878_variants.vcf.gz.tbi` files. After this, you can evaluate variants against the know truth set. We use the variants from GIAB consortium that can be found here: 

https://dovetail-public-data.s3-us-west-2.amazonaws.com/NA12878_data/NA12878_GIAB.vcf.gz

https://dovetail-public-data.s3-us-west-2.amazonaws.com/NA12878_data/NA12878_GIAB.vcf.gz.tbi

We also use a set of high-confident regions for evaluation. The bed file for evaluation can be found here:

https://dovetail-public-data.s3-us-west-2.amazonaws.com/NA12878_data/NA12878_highconf.bed


Once you download all the files, you can run the evaluation script as follows:

```
./snp_eval.sh NA12878_GIAB.vcf.gz NA12878_variants.vcf.gz NA12878_highconf.bed NA12878_variants_eval.txt
```

We can further phase these variants with HapCUT2 as follows:

```
./phasing_workflow.sh NA12878_variants.vcf.gz NA12878_OmniC.bam NA12878
```

The output will be stored in a file called `NA12878_hapcut_output.phased.VCF.gz`. The accuracy of phasing can be estimated if you have the truthset. Here we use the same GIAB truthset for evaluation:

```
./phasing_eval.sh NA12878_hapcut_output.phased.VCF.gz NA12878_GIAB.vcf.gz NA12878
```
