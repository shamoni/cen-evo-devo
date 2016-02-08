#! /usr/bin/env python
import Bio, argparse
from Bio import SeqIO

def main():
	parser = argparse.ArgumentParser(description='This script grabs fasta records whose ids match a specifier')
	parser.add_argument('input', action="store", type=str, help="Filename of a fasta file in current folder")
	parser.add_argument('-i','--id', action="store", type=str, help="Enter specifier to match ids")
	args = parser.parse_args()

	records = []
	for seq_record in SeqIO.parse(args.input, "fasta"):
		if args.id in seq_record.id:
			records.append(seq_record)

	output_handle = open(args.id+'_'+args.input, "w")
	SeqIO.write(records, output_handle, "fasta")

if __name__ == '__main__':
	main()
