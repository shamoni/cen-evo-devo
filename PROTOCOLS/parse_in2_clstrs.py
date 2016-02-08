#! /usr/bin/env python
import argparse, sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
	parser = argparse.ArgumentParser(description='This script splits a sequences in a fasta file by their cluster ids')
	parser.add_argument('input', action="store", type=str, help="Filename of a fasta file in current folder")
	args = parser.parse_args()

	records = []
	for record in SeqIO.parse(args.input, "fasta"):
		records.append(SeqRecord(record.seq.ungap("-"), id = record.id, description = ""))
	print >>sys.stderr, "Recovered %d sequences" %(len(records))

	clusters = set()
	for record in records:
		clusters.add(record.id[-1])

	record_by_cluster = [[] for i in range(len(clusters))]
	for record in records:
		clstr = int(record.id[-1])-1
		record_by_cluster[clstr].append(record)

	outdir = args.input.split(".")[0]+"_Clusters"
	os.mkdir(outdir)
	os.chdir(outdir)
	for i in range(len(clusters)):
		filename = args.input+"_Cluster"+str(i+1)+".fa"
		SeqIO.write(record_by_cluster[i],filename,"fasta")

if __name__ == '__main__':
	main()