#!/bin/bash

ref=$1
truth=$2
comp=$3
out=$4

SRCDIR=`dirname $0`

java -jar ${SRCDIR}/GenomeAnalysisTK.jar  -T GenotypeConcordance  -eval ${comp} -truth ${truth} -S ${out}
