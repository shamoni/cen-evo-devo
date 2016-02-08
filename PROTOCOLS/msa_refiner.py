#! /usr/bin/env python

import sys, argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

def del_all_gaps(alignment):
	rows = len(alignment)
	n = len(alignment[0])
	n_at_start = n
	i = 0
	while i < n:
		if alignment[:,i].count("-") == rows:
			if i == 0:
				alignment = alignment[:,1:]
			elif i+1 == n:
				alignment = alignment[:,:i]
			else:
				alignment = alignment[:,:i] + alignment[:,i+1:]
			n -= 1
		else:
			i += 1
	return alignment

def del_rows(alignment,t):
	rows = len(alignment)
	pos = len(alignment[0])
	gaps = []
	for i in range(pos):
		if rows - alignment[:,i].count("-") < t:
			gaps.append(i)

	to_remove = set()
	for c in gaps:
		for r in range(rows):
			if alignment[r,c] != "-":
				to_remove.add(alignment[r].id)

	filtered_records = []
	for record in alignment:
		if record.id not in to_remove:
			filtered_records.append(SeqRecord(record.seq, record.id, description=""))
	
	filtered_align = MultipleSeqAlignment(filtered_records)
	return filtered_align


def main():
	parser = argparse.ArgumentParser(description='This script adjusts a multiple sequence alignment \
		by removing sequences that insert rare gaps')
	parser.add_argument('input', action="store", type=str, help="Filename of .aln file in current folder")
	parser.add_argument('-t', '--threshold', action="store", type=int, default=20, help="Threshold to call rare insertions")
	args = parser.parse_args()

	input_file = args.input
	f = open(input_file)
	alignment = AlignIO.read(f, "fasta")
	print >>sys.stderr, "Recovered msa of lengh %d with %d aligned sequences" %(len(alignment[0]),len(alignment))

	t = args.threshold
	alignment = del_rows(alignment,t)
	print >>sys.stderr, "Kept %d aligned sequences after deleting sequences with insertions shared by < %d other sequences" %(len(alignment),t)

	alignment = del_all_gaps(alignment)
	print >>sys.stderr, "Trimmed msa has length %d with %d aligned sequences" %(len(alignment[0]),len(alignment))

	filename = input_file.split(".")[0]
	AlignIO.write(alignment, filename+".aln2", "fasta")

if __name__ == '__main__':
		main()

