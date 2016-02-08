#! /usr/bin/env python

import os, sys, math

#Generate a list with .fq files in current folder
ls = os.listdir(os.getcwd())
fqs = filter(lambda x: "R1.fq" in x, ls)

os.system("mkdir ../Merged_reads")
out = open('SeqPrep_results.txt', 'w')
out.write('File\tPairs Processed\tPairs Merged\tPairs With Adapters\tPairs Discarded\n')

for file in fqs:
        name = "_".join(file.split('_')[:-1])
        os.system("~/bin/SeqPrep-1.1/SeqPrep -f "+name+"_R1.fq \
                -r "+name+"_R2.fq -1"+name+"_F -2"+name+"_R \
                -q 30 -L 35 \
                -s"+name+"_merged.fq.gz 2>"+name+"_stdout.txt")
        #extract information from the _stdout file
        f = open(name+"_stdout.txt")
        lines = f.read().splitlines()[1:-1]
        proc = lines[0].split('\t')[1]
        merge = lines[1].split('\t')[1]
        adapt = lines[2].split('\t')[1]
        disc = lines[3].split('\t')[1]
        out.write(name+'\t'+proc+'\t'+merge+'\t'+adapt+'\t'+disc+'\n')

os.system("rm *_F")
os.system("rm *_R")
os.system("rm *_stdout.txt")
os.system("mv *_merged.fq.gz ../Merged_reads/")
os.system("mv SeqPrep_results.txt ../Merged_reads/")
