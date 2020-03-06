import pysam
import sys

bamfile = pysam.AlignmentFile(sys.argv[1],"rb")

stored_strand = {}
prev_chr = ''
prev_mapq = ''
line = ''
alignments = {}

chr_list = []

with open(sys.argv[2],'r') as f:
    for line in f:
        attrs = line.split()
        chr_list.append(attrs[0])

for r in bamfile.fetch(until_eof=True):
    #line += r.query_name+"\t
        mate_mapq = int(r.get_tag('MQ'))
        read_mapq = r.mapping_quality
        if mate_mapq < 60 or read_mapq < 60:                                        
                continue
        if r.is_supplementary or r.is_secondary:
                continue    
        if r.is_paired:
            if bamfile.get_reference_name(r.reference_id) in chr_list and bamfile.get_reference_name(r.next_reference_id) in chr_list:
                if chr_list.index(bamfile.get_reference_name(r.reference_id)) <= chr_list.index(bamfile.get_reference_name(r.next_reference_id)):
                    line = "1\t"
                    line += bamfile.get_reference_name(r.reference_id)+'\t'
                    line += str(r.reference_start)+'\t1\t'
                    line += '0\t'
                    line += bamfile.get_reference_name(r.next_reference_id)+'\t'
                    line += str(r.next_reference_start)+'\t0'   
                    print(line)
        line = ""
