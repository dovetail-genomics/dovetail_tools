#!/bin/bash

truth=$1
comp=$2
regions=$3
out=$4

SRCDIR=`dirname $0`

gatk Concordance  -eval ${comp} -truth ${truth} -L ${regions} -S ${out}
