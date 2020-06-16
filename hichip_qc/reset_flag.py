#!/usr/bin/env python3
import pysam
import sys

inputbam = sys.argv[1]

alignmentfile = pysam.AlignmentFile(sys.argv[1],'rb')
outbam  = pysam.AlignmentFile(sys.argv[2],'wb',template=alignmentfile)
for r in alignmentfile:
    if r.is_reverse:
        r.flag = 16
    else:
        r.flag = 0
    outbam.write(r)
