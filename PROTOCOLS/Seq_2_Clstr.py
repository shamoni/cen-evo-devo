#! /usr/bin/env python
import argparse, sys, os, csv
from Bio import SeqIO


def make_kmers(seq,ksize):
	x = []
	for i in range(len(seq) - ksize + 1):
		x.append(seq[i:i+ksize])
	return x

def slurp_kmers(filepath,ksize,cutoff):
	counts = {}
	pos = {}
	fp = open(filepath, "rU")
	for record in SeqIO.parse(fp, "fasta"):
		sequence = str(record.seq)
		sequence = sequence.upper()

		kmers = make_kmers(sequence,ksize)

		for i,kmer in enumerate(kmers):
			counts[kmer] = counts.get(kmer, 0) + 1
			if pos.has_key(kmer):
				pos[kmer].append(i)
			else:
				pos[kmer] = [i,]

	fp.close()
	keylist = counts.keys()
	for kmer in keylist:
		if counts[kmer] < cutoff:
			del counts[kmer]
			del pos[kmer]
	all_kmers = set(counts.keys())
	pos.update((key,list(set(value))) for key,value in pos.items())
	return all_kmers,pos

def parse_signatures(kmers):
	l = len(kmers)
	sigs = []
	for i in range(0,l):
		a = set(kmers[i])
		for j in range(0,l):
			if i != j:
				a = a - kmers[j]
			else:
				pass
		sigs.append(a)
	return sigs

def sigkmers(path,ksize,cutoff):
	ls = os.listdir(path)
	fasta = filter(lambda x: ".fa" in x,ls)
	fasta.sort()
	print >>sys.stderr, '\nProcessing Kmers from %d Clusters' %len(fasta)
	print >>sys.stderr, 'Identifying Kmers of size %d shared by at least %d sequences within the cluster' %(ksize,cutoff)
	print >>sys.stderr, ''

	k_by_clstr = []
	kpos_by_clstr = []
	for f in fasta:
		filepath = path+f
		kmers,pos = slurp_kmers(filepath,ksize,cutoff)
		k_by_clstr.append(kmers)
		kpos_by_clstr.append(pos)

	sig_k = parse_signatures(k_by_clstr)
	for i in range(0,len(fasta)): print >>sys.stderr, '%s; No. of Kmers %d; No. of Signature Kmers %d' \
												%(fasta[i],len(k_by_clstr[i]),len(sig_k[i]))
	print >>sys.stderr, ''

	sigkmer_positions = [{} for i in range(len(sig_k))]
	for i in range(len(sig_k)):
		for kmer in sig_k[i]:
			sigkmer_positions[i][kmer] = kpos_by_clstr[i][kmer]
	return (fasta,sig_k,sigkmer_positions)

def assign_read_to_cluster(read,read_rc,sig_kmers,ksize):
	x = set(make_kmers(read, ksize))
	y = set(make_kmers(read_rc, ksize))
	union = x.union(y)

	isect_list = []
	for index in range(len(sig_kmers)):
		cluster_specific_kmers = sig_kmers[index]
		if union.intersection(cluster_specific_kmers):
			isect_list.append(index)

	return isect_list

def assign_reads_to_clusters(filepath,sig_kmers,ksize,only_reads):
	n_zero = 0
	n_unique = 0
	n_multi = 0

	cluster_counts = [0]*len(sig_kmers)
	cluster_multi_counts = [0]*len(sig_kmers)
	unique_records = [[] for i in range(len(sig_kmers))]
	for n, record in enumerate(SeqIO.parse(open(filepath),"fasta")):
		if only_reads is not None and n > only_reads:
			break
		read = str(record.seq).upper()
		read_rc = str(record.seq.reverse_complement()).upper()
		isect_list = assign_read_to_cluster(read,read_rc,sig_kmers,ksize)

		if len(isect_list) == 0:
			n_zero += 1
		elif len(isect_list) == 1:
			n_unique += 1
			cluster_counts[isect_list[0]] += 1
			unique_records[isect_list[0]].append(record)
		else:
			assert len(isect_list) > 1
			n_multi += 1
			for i in isect_list:
				cluster_multi_counts[i] += 1

	summary = [n_zero,n_unique,n_multi]
	return summary, cluster_counts, cluster_multi_counts,unique_records

def get_cluster_distribution(path,sig_kmers,ksize,num_reads,outdir):
	ls = os.listdir(path)
	fasta = filter(lambda x: ".fa" in x,ls)
	fasta.sort()
	print >>sys.stderr, 'Processing reads from %d ChIPseq datasets' %len(fasta)
	
	data = {}
	for f in fasta:
		filepath = path+f
		data[f] = assign_reads_to_clusters(filepath,sig_kmers,ksize,num_reads)
		print >>sys.stderr, "Finished processing the file %s" %f

	return data

def write_data(data,outdir,clstr_id):
	os.mkdir(outdir)
	os.chdir(outdir)
	with open(outdir+"_summary.csv","wb") as outcsv:
		writer = csv.writer(outcsv)
		writer.writerow(["filename", "NOT_CEN","Unique_CEN","Multi_CEN"])
		for key,value in data.items():
			writer.writerow([key]+map(str,value[0]))

	with open(outdir+"_results.csv","wb") as outcsv:
		writer = csv.writer(outcsv)
		writer.writerow(["filename","Class"]+clstr_id)
		for key,value in data.items():
			writer.writerow([key]+["Unique"]+map(str,value[1]))
			writer.writerow([key]+["Multi"]+map(str,value[2]))

	for key,value in data.items():
		name = key.split('.')[0]
		for i in range(len(value[3])):
			clstr = clstr_id[i].split('.')[0]
			filename = name+"_"+clstr+".fa"
			SeqIO.write(value[3][i],filename,"fasta")

def write_signatures(sigkmer_positions,outdir,clstr_ids,k_len):
	with open(outdir+"_signatures.csv","wb") as outcsv:
		writer = csv.writer(outcsv)
		for i in range(len(clstr_ids)):
			for key,value in sigkmer_positions[i].items():
				writer.writerow([clstr_ids[i]]+[k_len]+[key]+map(str,value))

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-Fc", "--folder-w-clusters", dest="clusterpath", required=True)
	parser.add_argument("-Fr", "--folder-w-reads", dest="chipseqpath", required=True)
	parser.add_argument("-O", "--outdir", type=str,default="READ-2-CLSTR",dest="outdir")
	parser.add_argument("-K","--kmer-length", type=int, default=25, dest="k_len")
	parser.add_argument("-Kt","--cutoff-kmers",type=int,default=100, dest="k_co")
	parser.add_argument("-N","--max-no-of-reads-to-analyze",type=int,default=None, dest="num_reads")
	args = parser.parse_args()

	clstr_ids,sig_kmers,sigkmer_positions = sigkmers(args.clusterpath,args.k_len,args.k_co)
	data = get_cluster_distribution(args.chipseqpath,sig_kmers,args.k_len,args.num_reads,args.outdir)
	write_data(data,args.outdir,clstr_ids)
	write_signatures(sigkmer_positions,args.outdir,clstr_ids,args.k_len)

if __name__ == '__main__':
	main()
