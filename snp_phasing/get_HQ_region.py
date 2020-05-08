#!/usr/bin/env python3
import sys,argparse,pysam,os

# Create a bed file of the high confidence regions of the genome given a 
# particular alignment bam.  Here that is defined as all the regions where
# there are no MQ=0 reads.  

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-bam", required=True,help="Input bam file")
    parser.add_argument("-qthresh", required=False,type=int,default=0,help="Mapping quality < qthresh regions output")
    parser.add_argument("-bedroot",required=True,help="Root name for several output files produced. ")
    parser.add_argument("-contigsizes",required=True,help="Sizes of contigs.")
    args = parser.parse_args()
    return args   


def main(argv):
    args = parseArgs(argv)
    
    bam = pysam.AlignmentFile(args.bam)
    ofile = open(f"{args.bedroot}.bed",'w')
    for r in bam:
        if ((r.is_read1) and
            (not r.is_unmapped) and
            (not r.mate_is_unmapped) and
            (not r.is_supplementary) and
            (not r.is_secondary) and
            (r.mapping_quality < args.qthresh)):
                outstr = f"{r.reference_name}\t{r.reference_start}\t{r.reference_end}\t{r.mapping_quality}"
                ofile.write(outstr+'\n')    
                    
    
if __name__ == "__main__":
    main(sys.argv)
    
