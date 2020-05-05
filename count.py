import argparse
import pysam

parser = argparse.ArgumentParser()

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

def get_snr(coverage_in, coverage_out):
	coverage_in_peaks = 0
	linecount = 0

	with open(coverage_in,'r') as f:
		for line in f:
			attrs = line.split()
			coverage_in_peaks += float(attrs[-1])
			linecount += 1

	avg_coverage_in_peaks = coverage_in_peaks/linecount
	avg_coverage_in_peaks_r = float(str(coverage_in_peaks/linecount)[:5])

	coverage_outside_peaks = 0
	linecount = 0

	with open(coverage_out,'r') as f:
		for line in f:
			attrs = line.split()
			coverage_outside_peaks += float(attrs[-1])
			linecount += 1


	avg_coverage_outside_peaks = coverage_outside_peaks/linecount
	avg_coverage_outside_peaks_r = float(str(coverage_outside_peaks/linecount)[:5])

	snr_ratio = round(avg_coverage_in_peaks/avg_coverage_outside_peaks,2)

	return snr_ratio



parser.add_argument('-b1',help="bedfile1")
parser.add_argument('-b2',help="bedfile2")
parser.add_argument('-b3',help="bedfile3")
parser.add_argument('-b4',help="bedfile4")
parser.add_argument('-b5',help="bedfile5")
parser.add_argument('-cp1',help="covrage in peaks")
parser.add_argument('-co1',help="coverage outside peaks")
parser.add_argument('-cp2',help="covrage in peaks")
parser.add_argument('-co2',help="coverage outside peaks")
parser.add_argument('-cp3',help="covrage in peaks")
parser.add_argument('-co3',help="coverage outside peaks")
parser.add_argument('-cp4',help="covrage in peaks")
parser.add_argument('-co4',help="coverage outside peaks")
parser.add_argument('-prefix',help="prefix")
parser.add_argument('-bam',help="bamfile")
args = parser.parse_args()

in_peaks = get_read_count(args.b1)
in_500_peaks = get_read_count(args.b2)
in_1000_peaks = get_read_count(args.b3)
in_2000_peaks = get_read_count(args.b4)
in_5000_peaks = get_read_count(args.b5)

bamfile = pysam.AlignmentFile(args.bam,'rb')

paired_reads = 0

for r in bamfile.fetch(until_eof=True):
	 if r.is_paired:
	 	paired_reads += 1

total_valid_pairs = paired_reads//2
in_peaks_p = round(in_peaks *100.0/total_valid_pairs,2)
in_500_peaks_p = round(in_500_peaks *100.0/total_valid_pairs,2)
in_1000_peaks_p = round(in_1000_peaks *100.0/total_valid_pairs,2)
in_2000_peaks_p = round(in_2000_peaks *100.0/total_valid_pairs,2)
in_5000_peaks_p = round(in_5000_peaks *100.0/total_valid_pairs,2)



snr_ratio = get_snr(args.cp1,args.co1)
snr_ratio_1000 = get_snr(args.cp2,args.co2)
snr_ratio_2000 = get_snr(args.cp2,args.co3)
snr_ratio_5000 = get_snr(args.cp4,args.co4)

print(f"Total read pairs:\t{total_valid_pairs}")
print(f"Total read pairs in peaks:\t{in_peaks}({in_peaks_p}%)")
print(f"Total read pairs in 500 bp around peaks:\t{in_500_peaks}({in_500_peaks_p}%)")
print(f"Total read pairs in 1000 bp around peaks:\t{in_1000_peaks}({in_1000_peaks_p}%)")
print(f"Total read pairs in 2000 bp around peaks:\t{in_2000_peaks}({in_2000_peaks_p}%)")
print(f"Total read pairs in 5000 bp around peaks:\t{in_5000_peaks}({in_5000_peaks_p}%)")
print(f"Signal to noise ratio for reads in peaks:\t{snr_ratio}")
print(f"Signal to noise ratio for reads in 1000 bp around peaks:\t{snr_ratio_1000}")
print(f"Signal to noise ratio for reads in 2000 bp around peaks:\t{snr_ratio_2000}")
print(f"Signal to noise ratio for reads in 5000 bp around peaks:\t{snr_ratio_5000}")
