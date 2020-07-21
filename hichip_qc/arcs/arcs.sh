#!/bin/bash


mcool=$1
bamfile=$2
outputprefix=$3

cooler dump --join ${mcool}::/resolutions/1000 > ${outputprefix}_HiChIP.cool.txt 
bamCoverage -b ${bamfile} -of bedgraph -p 36 -o ${outputprefix}.coverage.bedgraph



