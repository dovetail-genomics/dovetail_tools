#!/usr/bin/env python3

import argparse

'''
Mapping Quality Threshold         : 40
Read1                             : 70,844,596
Read2                             : 70,752,654
Mapped pairs                      : 69,493,638
PCR dupe pairs                    : 8,693,432
Mapped nondupe pairs              : 60,800,206
Valid Pairs (cis>1000bp + trans)  : 30,171,190
Mapped nondupe pairs cis          : 50,249,773
Mapped nondupe pairs cis <=1000bp : 30,629,016
Mapped nondupe pairs cis >1000bp  : 19,620,757
Mapped nondupe pairs cis >10000bp : 6,227,109
Mapped nondupe trans pairs        : 10,550,433
Expected unique pairs at 300M sequencing:  218463642.2
Total ChIP peaks:	286,152
Mean ChIP peak size:	198 bp
Median ChIP peak size:	204 bp
Total read pairs in peaks:	1,369,652(1.88%)
Total read pairs in 500 bp around peaks:	10,654,570(14.6%)
Total read pairs in 1000 bp around peaks:	18,747,602(25.69%)
Total read pairs in 2000 bp around peaks:	31,348,497(42.95%)
Total read pairs in 5000 bp around peaks:	51,324,571(70.32%)
'''

parser = argparse.ArgumentParser()
parser.add_argument('-i', help="Input file")
parser.add_argument('-r', help="r1 fastq")
parser.add_argument("-o", help="Output file")
args = parser.parse_args()

header = ["Library_Name"]
body = [args.r]
ofile = open(args.o,'w')

with open(args.i, 'r') as f:
    for line in f:
        attrs = line.strip().split(":")
        header.append(attrs[0].strip())
        body.append(attrs[1].strip())

header_str = "\t".join(header)
body_str = "\t".join(body)

ofile.write(header_str+'\n'+body_str)
ofile.close()
