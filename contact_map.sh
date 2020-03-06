#!/bin/bash

bam=$1
outputfile=$2
asm=$3
tmpdir=/tmp
prefix=$4
#prefix=`basename  ${bam} | sed 's/.bam//g'`
SCRIPT_PATH=`realpath $0`
SCRIPT_PATH=`dirname $SCRIPT_PATH`

samtools faidx $asm
cut -f 1,2 $asm.fai > /tmp/chr_size.txt
chrlength="/tmp/chr_size.txt"
python ${SCRIPT_PATH}/convert_alignments_to_tsv.py ${bam} ${chrlength} \
    | awk '{if ($2 > $6) {print $1"\t"$6"\t"$7"\t"$8"\t"$5"\t"$2"\t"$3"\t"$4} else {print}}' \
    | sort -k2,2d -k6,6d -T ${tmpdir} --parallel=40 \
    > ${tmpdir}/${prefix}_juicer_alignments.txt

java -jar -Xmx30g \
     ${SCRIPT_PATH}/juicertools.jar pre \
    ${tmpdir}/${prefix}_juicer_alignments.txt ${outputfile} ${chrlength}

hic2cool convert -r 1000  ${outputfile} ${prefix}.cool

cooler zoomify --balance -p 40 ${prefix}.cool
