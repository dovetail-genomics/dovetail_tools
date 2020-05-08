#!/bin/bash

hapcut_phasing=$1
truth_phasing=$2
out_prefix=$3

whatshap compare --only-snvs --switch-error-bed ${out_prefix}_switch_erros.txt --tsv-pairwise ${out_prefix}_phasing_eval.txt ${truth_phasing} ${hapcut_phasing} 
