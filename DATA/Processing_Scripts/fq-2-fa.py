#! /usr/bin/env python

import os

#Generate a list with .fq files in current folder
ls = os.listdir(os.getcwd())
fqs = filter(lambda x: ".fq" in x and ".clstr" not in x, ls)

for fq_file in fqs:
	f = open(fq_file, "r")
	name = fq_file.split(".")[0]
	o = open(name+".fa", "w")
	while True:
        	seqid = f.readline()
        	seqid = seqid.split(' ')[0]
        	if seqid == "":
            		break
         	seq = f.readline()
         	plus = f.readline()
         	qual = f.readline()
		o.write(">"+seqid[1:]+"\n"+seq)
	f.close()
	o.close()
