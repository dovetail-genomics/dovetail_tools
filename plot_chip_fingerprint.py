#!/usr/bin/env python3

import matplotlib.pyplot as plt
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-table", help="count table")
parser.add_argument("-output", help="output plot")

args = parser.parse_args()

f = open(args.table,'r')
lines = f.readlines()
sample = lines[1]
counts = []
fraction_wrt_max = []
for line in lines[2:]:
    counts.append(int(line.strip()))

max_coverage = np.max(counts)

paired_list = []
for count in counts:
    fraction_wrt_max = count/max_coverage
    paired_list.append([fraction_wrt_max, count])


sorted_pairs_list = sorted(paired_list, key = lambda x:x[1])
csum = 0

x = []
y = []
coverage_sum = np.sum(counts)
for i in range(1, len(sorted_pairs_list)):
    sorted_pairs_list[i][1] += sorted_pairs_list[i-1][1]
    x.append(sorted_pairs_list[i][1]/coverage_sum)
    y.append(i/len(sorted_pairs_list))

plt.plot(y,x)
plt.grid()
plt.title("ChIP fingerprint")
plt.xlabel("Fraction of bins")
plt.ylabel("Coverage wrt highest coverage bin")
plt.savefig(args.output,dpi=300)
