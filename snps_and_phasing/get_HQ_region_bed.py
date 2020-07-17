#!/usr/bin/env python3
import os,sys,argparse,pysam

# Normally bwa uses the tight read pair gap size distribution to place some 
# otherwise ambiguously mapped mates of mapped reads.  By nature, 
# OmniC has a very broad read pair gap size distribution which differs 
# from the expectation of bwa.  As a result, some fraction of mate pairs 
# can not be placed unambiguously and this uncertainty hurts SNP calling in 
# some difficult to map regions.  It might be possible to alter bwa to better
# fit the expectations of OmniC data, but in the mean time we can 
# increase our confidence in called SNPs by restricting ourselves to 
# regions where no reads have MQ0.  This is a stringent criteria and 
# can be loosened in various ways if needed.  

# This script creates a bed file that describes all such regions.  This 
# bed file can then be used to filter a vcf to produce the high-confidence
# vcf. 

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-bam", required=True,help="Input bam file (indexed sorted)")
    parser.add_argument("-bedroot",required=True,help="Root name for several output files. ")
    parser.add_argument("-qthresh", required=False,type=int,default=0,help="Mapping quality < qthresh regions output")
    parser.add_argument("-region", required=False, default=None)
    args = parser.parse_args()
    return args   
            
def main(argv):
    args = parseArgs(argv)
    
    bam = pysam.AlignmentFile(args.bam)
    with open(f"{args.bedroot}_MQ0.bed","w") as bout:
        for r in bam.fetch(region=args.region):
            if ((r.is_read1) and
                (not r.is_unmapped) and
                (not r.mate_is_unmapped) and
                (not r.is_supplementary) and
                (not r.is_secondary) and
                (r.mapping_quality <= args.qthresh)):
                    outstr = f"{r.reference_name}\t{r.reference_start}\t{r.reference_end}\t{r.mapping_quality}"
                    bout.write(outstr+"\n")    
            
    cmd=f"bedtools merge -i {args.bedroot}_MQ0.bed > {args.bedroot}_merged.bed"
    print(cmd)
    os.system(cmd)
    
    # Determine genome contig sizes from the alignment
    sizesfile = f"{args.bam}.stats.txt"
    cmd = f"samtools idxstats {args.bam} | grep -v '*' > {sizesfile}"
    print(cmd)
    os.system(cmd)
    
    # Complement LQ bed relative to the genome contigs to get a high confidence bed. 
    cmd=f"bedtools complement -i {args.bedroot}_merged.bed -g {sizesfile} > {args.bedroot}_highconf.bed"
    print(cmd)
    os.system(cmd)  

    # Remove MQ0 file as it's large and unnecessary after merge. 
    cmd = f"rm -f {args.bedroot}_MQ0.bed"
    print(cmd)
    os.system(cmd)    
    
if __name__ == "__main__":
    main(sys.argv)
