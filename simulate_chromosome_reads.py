#! /bin/python3

from Bio import SeqIO
import math
import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--ref')
parser.add_argument('--k', type=int)
parser.add_argument('--chr')
args = parser.parse_args()

chr_index = SeqIO.index(args.ref, 'fasta')

def list_kmers(chr_sequence, k):
    kmer_list = []
    if len(chr_sequence) - k + 1 > 0:
        for start in range(len(chr_sequence))[: len(chr_sequence) - k + 1]:
            current_kmer = chr_sequence[start:start+k]
            if current_kmer != 'N'*k:        
                kmer_list.append(current_kmer)
    return kmer_list

chr_sequence = str(chr_index[args.chr].seq)

kmer_list = list_kmers(chr_sequence, args.k)
digits = int(math.log10(len(kmer_list))) + 1

filename='header.sam'
with open(filename, 'w') as prebam:
    prebam.write('%s\t%s\n%s\t%s\t%s\n' % ('@HD','VN:1.9', '@RG','ID:simulated','SM:simulated'))

counter = 0
filename = args.chr + '.ubam'
header_file = pysam.AlignmentFile('header.sam', check_sq = False)
with pysam.AlignmentFile(filename, 'wb', template=header_file) as outf:
    for k in kmer_list[:1000]:
        seq = pysam.AlignedSegment()
        seq.query_name = str(counter).zfill(digits)
        seq.query_sequence=k
        seq.flag = 4
        seq.query_qualities =  pysam.qualitystring_to_array(''.join('*' * args.k))
        outf.write(seq)
        counter = counter+1
