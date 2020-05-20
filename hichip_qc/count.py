#!/usr/bin/env python3

import argparse
import pysam
import numpy as np

def get_read_count(bedfile):
    count = 0
    last_read = ''
    with open(bedfile,'r') as f:
        for line in f:
            attrs = line.split()
            read = attrs[4].split('/')[0]
            if read != last_read:
                count += 1
            last_read = read
    return count


def get_snr(in_peaks, paired_reads, bp_in_loops , ref_length, number_of_loops, padding):
    loop_in_length = 2*padding*number_of_loops + bp_in_loops
    outside_loop_length = ref_length - loop_in_length
    coverage_in_loop = 2*150*in_peaks/loop_in_length
    coverage_outside_loop = 2*150*(paired_reads - in_peaks)/outside_loop_length
    snr_ratio = round(coverage_in_loop/coverage_outside_loop,2)
    return snr_ratio


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b1',help="bedfile1")
    parser.add_argument('-b2',help="bedfile2")
    parser.add_argument('-b3',help="bedfile3")
    parser.add_argument('-b4',help="bedfile4")
    parser.add_argument('-b5',help="bedfile5")
    parser.add_argument('-peaks',help="input peaks")
    parser.add_argument('-bam',help="bamfile")
    args = parser.parse_args()

    in_peaks = get_read_count(args.b1)
    in_500_peaks = get_read_count(args.b2)
    in_1000_peaks = get_read_count(args.b3)
    in_2000_peaks = get_read_count(args.b4)
    in_5000_peaks = get_read_count(args.b5)

    in_500_peaks_fmt = format(in_500_peaks,",d")
    in_1000_peaks_fmt = format(in_1000_peaks,",d")
    in_2000_peaks_fmt = format(in_2000_peaks,",d")
    in_5000_peaks_fmt = format(in_5000_peaks, ",d")

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

    total_valid_pairs = paired_reads//2
    in_peaks_p = round(in_peaks *100.0/total_valid_pairs,2)
    in_500_peaks_p = round(in_500_peaks *100.0/total_valid_pairs,2)
    in_1000_peaks_p = round(in_1000_peaks *100.0/total_valid_pairs,2)
    in_2000_peaks_p = round(in_2000_peaks *100.0/total_valid_pairs,2)
    in_5000_peaks_p = round(in_5000_peaks *100.0/total_valid_pairs,2)
    median_peak_size = format(int(np.median(peak_size)),",d")
    mean_peak_size = format(int(np.mean(peak_size)), ",d")

    
    print(f"Total ChIP peaks:\t{number_of_loops}")
    print(f"Mean ChIP peak size:\t{mean_peak_size} bp")
    print(f"Median ChIP peak size:\t{median_peak_size} bp")
    print(f"Total read pairs in 500 bp around peaks:\t{in_500_peak_fmt}({in_500_peaks_p}%)")
    print(f"Total read pairs in 1000 bp around peaks:\t{in_1000_peaks_fmt}({in_1000_peaks_p}%)")
    print(f"Total read pairs in 2000 bp around peaks:\t{in_2000_peaks_fmt}({in_2000_peaks_p}%)")
    print(f"Total read pairs in 5000 bp around peaks:\t{in_5000_peaks_fmt}({in_5000_peaks_p}%)")
