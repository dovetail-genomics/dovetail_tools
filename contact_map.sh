#!/bin/bash

pairs=$1
ref=$2
prefix=$3
cores=$4
#prefix=`basename  ${bam} | sed 's/.bam//g'`
SCRIPT_PATH=`realpath $0`
SCRIPT_PATH=`dirname $SCRIPT_PATH`

#create index of reference and chromosome size file
samtools faidx $ref
cut -f 1,2 $ref.fai > ${prefix}.genome
chrlength="${prefix}.genome"

#convert the pairs  file into juicer compatible text file.
zcat ${pairs} \
	| grep -v "#" \
	| grep "UU" \
	| awk '{print "1\t"$2"\t"$3"\t1\t0\t"$4"\t"$5"\t0"}' \
	| awk '{if ($2 > $6) {print $1"\t"$6"\t"$7"\t"$8"\t"$5"\t"$2"\t"$3"\t"$4} else {print}}' \
	| sort -k2,2d -k6,6d  --parallel=${cores} > ${prefix}_juicer_alignments.txt

#Convert the sorted text file into .hic file
java -jar -Xmx30g \
     ${SCRIPT_PATH}/juicertools.jar pre \
    ${prefix}_juicer_alignments.txt ${prefix}.hic ${chrlength}


#now convert pairs to .hic
pairix ${pairs}

cooler cload pairix \
	 -p ${cores} \
	 ${prefix}.genome:1000 \
	 ${pairs} \
	 ${prefix}.cool

#generate multi-resolution cooler file from 1kb resolution cooler file
cooler zoomify --balance -p ${cores}  ${prefix}.cool
