#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
from tabulate import tabulate

def get_read_count(bedfile):
    count = 0
    last_read = ''
    prev_id = ''
    reads = set()
    with open(bedfile,'r') as f:
        for line in f:
            attrs = line.split()
            read = attrs[3]
            reads.add(read)
    return len(reads)


def get_expected_observed(num_peaks, window_size, ref_length, reads_total, reads_observed):
    prob_read = num_peaks*window_size/ref_length
    expected_reads = prob_read*reads_total
    observed_expected_ratio = reads_observed/expected_reads
    return observed_expected_ratio


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b1',help="bedfile1")
    parser.add_argument('-b2',help="bedfile2")
    parser.add_argument('-b3',help="bedfile3")
    parser.add_argument('-b4',help="bedfile4")
    parser.add_argument('-b5', help='bedfile5')
    parser.add_argument('-peaks',help="input peaks")
    parser.add_argument('-bam',help="bamfile")
    args = parser.parse_args()
    

    in_peaks = get_read_count(args.b1)
    in_500_peaks = get_read_count(args.b2)
    in_1000_peaks = get_read_count(args.b3)
    in_2000_peaks = get_read_count(args.b4)
    in_blacklist = get_read_count(args.b5)

    in_peaks_fmt = format(in_peaks,",d")
    in_500_peaks_fmt = format(in_500_peaks,",d")
    in_1000_peaks_fmt = format(in_1000_peaks,",d")
    in_2000_peaks_fmt = format(in_2000_peaks,",d")
    in_blacklist_fmt = format(in_blacklist,",d")

    bamfile = pysam.AlignmentFile(args.bam,'rb')

    paired_reads = 0


    bp_in_loops = 0

    number_of_loops = 0
    peak_size = []
    with open(args.peaks,'r') as f:
        for line in f:
            attrs = line.strip().split('\t')
            bp_in_loops += int(attrs[2]) - int(attrs[1])
            peak_size.append(int(attrs[2]) - int(attrs[1]))
            number_of_loops += 1

    ref_length = sum(bamfile.lengths)

    bp_outside_loops = ref_length - bp_in_loops
    
    total_no_reads = 0
    for r in bamfile.fetch(until_eof=True):
         if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
            total_no_reads += 1

    
    ratio_in_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size), ref_length, total_no_reads, in_peaks), 2)
    ratio_in_500_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size) + 500, ref_length, total_no_reads, in_500_peaks),2)
    ratio_in_1000_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size)+ 1000, ref_length, total_no_reads, in_1000_peaks),2)
    ratio_in_2000_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size)+ 2000, ref_length, total_no_reads, in_2000_peaks),2)
    

    in_peaks_p = round(in_peaks *100.0/total_no_reads,2)
    in_500_peaks_p = round(in_500_peaks *100.0/total_no_reads,2)
    in_1000_peaks_p = round(in_1000_peaks *100.0/total_no_reads,2)
    in_2000_peaks_p = round(in_2000_peaks *100.0/total_no_reads,2)
    in_blacklist_p = round(in_blacklist*100.0/total_no_reads,2)

    median_peak_size = format(int(np.median(peak_size)),",d")
    mean_peak_size = format(int(np.mean(peak_size)), ",d")
    number_of_loops = format(number_of_loops, ",d")
   
    table = []
    table.append(["Total ChIP peaks", number_of_loops])
    table.append(["Mean ChIP peak size", f"{mean_peak_size} bp"])
    table.append(["Median ChIP peak size", f"{mean_peak_size} bp"])
    table.append(["Total reads in blacklist regions", f"{in_blacklist_fmt}", f"{in_blacklist_p}%"])
    table.append(["Total reads  in peaks", f"{in_peaks_fmt}", f"{in_peaks_p}%"])
    table.append(["Total reads in 500 bp around peaks", in_500_peaks_fmt, f"{in_500_peaks_p}%"])
    table.append(["Total reads in 1000 bp around peaks", in_1000_peaks_fmt, f"{in_1000_peaks_p}%"])
    table.append(["Total reads in 2000 bp around peaks", in_2000_peaks_fmt, f"{in_2000_peaks_p}%"])
    table.append(["Observed/Expected ratio for reads in peaks", ratio_in_peaks])
    table.append(["Observed/Expected ratio for reads in 500bp around", ratio_in_500_peaks])
    table.append(["Observed/Expected ratio for reads in 1000bp around  peaks", ratio_in_1000_peaks])
    table.append(["Observed/Expected ratio for reads in 2000bp around  peaks", ratio_in_2000_peaks])
    print(tabulate(table,tablefmt="plain"))
