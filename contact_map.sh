#!/bin/bash

bam=$1
ref=$2
tmpdir=/tmp
prefix=$3
#prefix=`basename  ${bam} | sed 's/.bam//g'`
SCRIPT_PATH=`realpath $0`
SCRIPT_PATH=`dirname $SCRIPT_PATH`

#create index of reference and chromosome size file
samtools faidx $ref
cut -f 1,2 $ref.fai > /tmp/chr_size.txt
chrlength="/tmp/chr_size.txt"

#convert the bam file into juicer compatible text file.
#Make sure it's sorted appropriately in the format juicer requires
python ${SCRIPT_PATH}/convert_alignments_to_tsv.py ${bam} ${chrlength} \
    | awk '{if ($2 > $6) {print $1"\t"$6"\t"$7"\t"$8"\t"$5"\t"$2"\t"$3"\t"$4} else {print}}' \
    | sort -k2,2d -k6,6d -T ${tmpdir} --parallel=40 \
    > ${tmpdir}/${prefix}_juicer_alignments.txt

#Convert the sorted text file into .hic file
java -jar -Xmx30g \
     ${SCRIPT_PATH}/juicertools.jar pre \
    ${tmpdir}/${prefix}_juicer_alignments.txt ${prefix}.hic ${chrlength}

#convert this generated hic file into cool file
hic2cool convert -r 1000  ${outputfile} ${prefix}.cool

#generate multi-resolution cooler file from 1kb resolution cooler file
cooler zoomify --balance -p 4 ${prefix}.cool
