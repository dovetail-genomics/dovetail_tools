#!/usr/bin/env python3

import argparse
import pysam
import numpy as np

def get_read_count(bedfile):
    count = 0
    last_read = ''
    prev_id = ''
    with open(bedfile,'r') as f:
        for line in f:
            attrs = line.split()
            read, read_id = attrs[3].split('/')
            if read == last_read and read_id != prev_id:
                count += 1
            last_read = read
            prev_id = read_id
    return count


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
    parser.add_argument('-peaks',help="input peaks")
    parser.add_argument('-bam',help="bamfile")
    args = parser.parse_args()
    

    in_peaks = get_read_count(args.b1)
    in_500_peaks = get_read_count(args.b2)
    in_1000_peaks = get_read_count(args.b3)
    in_2000_peaks = get_read_count(args.b4)

    in_peaks_fmt = format(in_peaks,",d")
    in_500_peaks_fmt = format(in_500_peaks,",d")
    in_1000_peaks_fmt = format(in_1000_peaks,",d")
    in_2000_peaks_fmt = format(in_2000_peaks,",d")

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

    for r in bamfile.fetch(until_eof=True):
         if r.is_paired:
            paired_reads += 1

    
    ratio_in_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size), ref_length, paired_reads, in_peaks), 2)
    ratio_in_500_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size) + 500, ref_length, paired_reads, in_500_peaks),2)
    ratio_in_1000_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size)+ 1000, ref_length, paired_reads, in_1000_peaks),2)
    ratio_in_2000_peaks = round(get_expected_observed(number_of_loops,np.mean(peak_size)+ 2000, ref_length, paired_reads, in_2000peaks),2)
    
    total_valid_pairs = paired_reads//2
    in_peaks_p = round(in_peaks *100.0/total_valid_pairs,2)
    in_500_peaks_p = round(in_500_peaks *100.0/total_valid_pairs,2)
    in_1000_peaks_p = round(in_1000_peaks *100.0/total_valid_pairs,2)
    in_2000_peaks_p = round(in_2000_peaks *100.0/total_valid_pairs,2)
    median_peak_size = format(int(np.median(peak_size)),",d")
    mean_peak_size = format(int(np.mean(peak_size)), ",d")
    number_of_loops = format(number_of_loops, ",d")
    
    print(f"Total ChIP peaks:\t{number_of_loops}")
    print(f"Mean ChIP peak size:\t{mean_peak_size} bp")
    print(f"Median ChIP peak size:\t{median_peak_size} bp")
    print(f"Total read pairs in peaks:\t{in_peaks_fmt}({in_peaks_p}%)")
    print(f"Total read pairs in 500 bp around peaks:\t{in_500_peaks_fmt}({in_500_peaks_p}%)")
    print(f"Total read pairs in 1000 bp around peaks:\t{in_1000_peaks_fmt}({in_1000_peaks_p}%)")
    print(f"Total read pairs in 2000 bp around peaks:\t{in_2000_peaks_fmt}({in_2000_peaks_p}%)")
    print(f"Observed/Expected ratio for reads in peaks:\t{ratio_in_peaks}")
    print(f"Observed/Expected ratio for reads in 500bp around  peaks:\t{ratio_in_500_peaks}")
    print(f"Observed/Expected ratio for reads in 1000bp around  peaks:\t{ratio_in_1000_peaks}")
    print(f"Observed/Expected ratio for reads in 2000bp around  peaks:\t{ratio_in_2000_peaks}")
