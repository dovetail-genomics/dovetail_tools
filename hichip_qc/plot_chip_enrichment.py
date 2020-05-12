#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-bam", help="Input BAM file")
parser.add_argument("-peaks", help="ChiSeq peaks in encode format")
parser.add_argument("-output", help="Output file")

args = parser.parse_args()

peak_data = pd.read_csv(args.peaks, sep="\t", header=None, keep_default_na=False)
peak_data.columns = ["chromosome", "start", "end", "A", "B", "C", "Signal_value", "D", "E", "offset"]
Q1 = peak_data['Signal_value'].quantile(0.25)
Q3 = peak_data['Signal_value'].quantile(0.75)
IQR = Q3 - Q1

peak_data_filtered = peak_data.query('(@Q3 + 1.5 * @IQR) <= Signal_value')


count = 0

coverage = dict()

for i in range(-10000, 10001):
    coverage[i] = 0

for index, row in peak_data_filtered.iterrows():
    chrom = row["chromosome"]
    center = row["start"] + row["offset"]
    start = center - 10000
    end = center + 10000
    x = []
    y = []
    index = -10000
    proc = subprocess.Popen(['samtools','mpileup', '-A', args.bam, '-r', f"{chrom}:{start}-{end}"],stdout=subprocess.PIPE)
    while True:
        line = proc.stdout.readline().decode('ascii')
        if not line: 
            break
        attrs = line.split('\t')
        coverage[index] += int(attrs[3])
        index += 1
    count += 1
    if count ==100:
        break
x = list(coverage.keys())
y = list(coverage.values())
N = 200
y_smooth = np.convolve(y, np.ones((N,))/N, mode='valid')
y_smooth_norm =[float(i)/np.mean(y_smooth) for i in y_smooth]
x_ticks = [0, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000] 
x_labels = ["-10kb", "-7.5kb", "-5kb", "-2.5kb", "0kb", "+2.5kb", "+5kb", "+7.5kb", "+10kb"]

plt.xticks(ticks=x_ticks, labels=x_labels)
plt.grid()
#plt.tight_layout()
plt.plot(y_smooth_norm)
plt.title("Coverage around ChIP peaks")
plt.xlabel("Distance from the peak center")
plt.ylabel("Fold coverage change  based on average coverage")
plt.savefig(args.output, dpi=300)
