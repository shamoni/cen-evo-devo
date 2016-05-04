#!/usr/bin/python

import os, sys, math, argparse

execfile('/software/modules/3.2.10/x86_64-linux-ubuntu14.04/Modules/3.2.10/init/python.py')
module('load','bwa')
module('load','samtools')

parser = argparse.ArgumentParser(description='Program to map reads using bwa')
parser.add_argument('-f', action="store", dest="pathFolder", required = True, help="Path to folder containing fastq flies")
parser.add_argument('-d', action="store", dest="pathGenome", required = True, help="Path to indexed genome with the index prefix")
parser.add_argument('-o',action="store",dest="outdir",required = True, help = "name of outdirectory to place all the files")
args = parser.parse_args()

pathFolder = args.pathFolder
pathGenome = args.pathGenome
outdir = args.outdir

print(pathGenome)

li = os.listdir(pathFolder)
fqs = filter(lambda x: ".fq" in x, li)
index = pathGenome.split('/')[-1]

os.system("mkdir "+outdir)

for file in fqs:
	name = "_".join(file.split('_')[:-1])
	outname = name+"_mapto_"+index
	os.system("bwa aln -t 10 "+pathGenome+" "+pathFolder+file+" > "+outname+".sai")
	os.system("bwa samse -n 10 -f "+outname+".sam "+pathGenome+" "+outname+".sai "+pathFolder+file)
	os.system("samtools view -b -S "+outname+".sam >"+outname+"_unsorted.bam")
	os.system("samtools sort "+outname+"_unsorted.bam "+outname)
	os.system("samtools index "+outname+".bam "+outname+".bai")
	os.system("rm "+outname+"_unsorted.bam")
	os.system("mv "+outname+"* "+outdir+"/")
