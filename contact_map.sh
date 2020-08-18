#!/bin/bash

if [  $# -le 3 ]
then
        echo "Too few arguments. Please provide all the required arguments."
        echo "Usage: ./contact_map.sh <pair.gz file> <reference_fasta> <output_prefix>> <num_cores>"
        exit 1
fi



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
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' ${pairs} -o ${prefix}.filtered.pairs.gz

zcat ${prefix}.filtered.pairs.gz \
	| grep -v "#" \
	| awk '{print "1\t"$2"\t"$3"\t1\t0\t"$4"\t"$5"\t0"}' \
	| awk '{if ($2 > $6) {print $1"\t"$6"\t"$7"\t"$8"\t"$5"\t"$2"\t"$3"\t"$4} else {print}}' \
	| sort -k2,2d -k6,6d  --parallel=${cores} > ${prefix}_juicer_alignments.txt

#Convert the sorted text file into .hic file
java -jar -Xmx30g \
     ${SCRIPT_PATH}/juicertools.jar pre \
    ${prefix}_juicer_alignments.txt ${prefix}.hic ${chrlength}


#now convert pairs to .hic
pairix ${prefix}.filtered.pairs.gz

cooler cload pairix \
	 -p ${cores} \
	 ${prefix}.genome:1000 \
	 ${prefix}.filtered.pairs.gz \
	 ${prefix}.cool

#generate multi-resolution cooler file from 1kb resolution cooler file
cooler zoomify --balance -p ${cores}  ${prefix}.cool

zcat ${prefix}.filtered.pairs.gz | \
	grep -v '#' | \
	awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | \
	gzip -c > ${prefix}.hicpro.valid.pairs.gz
