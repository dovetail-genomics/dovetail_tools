#!/usr/bin/env python3
import sys,os,csv,gzip

# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  omnic37X
# chr20   61791   .       C       A       232.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=-1.577e+00;DP=169;ExcessHet=3.0103;FS=55.090;MLEAC=1;MLEAF=0.500;MQ=52.91;MQRankSum=-3.104e+00;QD=4.23;ReadPosRankSum=2.82;SOR=0.185    GT:AD:DP:GQ:PGT:PID:PL:PS 0|1:46,9:169:99:0|1:61791_C_A:240,0,1730:61791

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-vcfin", required=True,help="VCF in")
    parser.add_argument("-qdthresh", required=False,type=float,default=2,help="QD threshold.  SNPS with QD<= qdthresh omitted.")
    parser.add_argument("-fsthresh", required=False,type=float,default=2,help="QD threshold.  SNPS with QD<= qdthresh omitted.")
    parser.add_argument("-vcfout",required=True,help="VCF output name (ending in .vcf)"
    args = parser.parse_args()
    return args   

def getInfoDict(line):
    """Get a dictionary of the value pairs in the info field"""
    infodict = dict()
    columns = line.split("\t")
    info = columns[7]
    for item in info.split(";"):
        fields = item.split("=")
        if len(fields) == 2:
            infodict[fields[0]] = fields[1]
                
    return infodict

def main(argv):
    args = parseArgs(argv)   
 
    with open(args.vcfout,'w') as fout:
        with gzip.open(args.vcfin,'rt') as fin:
            for line in fin:
                if line.startswith("#"):
                    fout.write(line)
                else:
                    infodict = getInfoDict(line)
                    if float(infodict['QD']) >= args.qdthresh and float(infodict['FS']) <= args.fsthresh:
                        fout.write(line)
    
    print(f"bgzip {args.vcfout}")
    os.system(f"bgzip {args.vcfout}")
    print(f"tabix -f -p vcf {args.vcfout}.gz")
    os.system(f"tabix -f -p vcf {args.vcfout}.gz")
    

if __name__ == "__main__":
    main(sys.argv)
