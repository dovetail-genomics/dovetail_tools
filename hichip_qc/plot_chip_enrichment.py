#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import subprocess

'''
chr1    826678  826954  Peak_12342      156     .       6.23321 12.01192        9.02289 152
chr1    958283  959249  Peak_589        393     .       10.52320        28.16420        23.14729        442
chr1    1115852 1116138 Peak_9793       186     .       6.88934 14.01684        10.81947        155
chr1    1157596 1158571 Peak_7928       209     .       6.83383 15.68351        12.29603        179
chr1    1230531 1231969 Peak_3912       272     .       7.32373 19.95536        16.04756        988
chr1    1324817 1326133 Peak_3564       280     .       8.85772 20.48661        16.51214        554
'''

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

for i in range(-1000, 1001):
    coverage[i] = 0

count = 0
for num, row in peak_data_filtered.iterrows():
    chrom = row["chromosome"]
    center = row["start"] + row["offset"]
    start  = center - 1000
    end = center + 1000
    x = []
    y = []
    count += 1
    proc = subprocess.Popen(['samtools','mpileup', '-A', args.bam, '-r', f"{chrom}:{start}-{end}"],stdout=subprocess.PIPE)
    while True:
        line = proc.stdout.readline().decode('ascii')
        if not line: 
            break
        attrs = line.split('\t')
        pos = int(attrs[1])
        coverage[pos - center] += int(attrs[3])

x = list(coverage.keys())
y = list(coverage.values())
x_ticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000] 
x_labels = ["-1000 bp", "-750 bp", "-500 bp", "-250 bp", "0 bp", "+250 bp", "+500 bp", "+750 bp", "+1000 bp"]

plt.xticks(ticks=x_ticks, labels=x_labels, fontsize=7)
plt.grid()
plt.plot(y)
plt.title("Coverage around ChIP peaks")
plt.xlabel("Distance from the peak center")
plt.ylabel("Fold coverage change  based on average coverage")
plt.savefig(args.output, dpi=300)
