#! /usr/bin/env python

import os, sys, math

#Generate a list with .fq files in current folder
ls = os.listdir(os.getcwd())
fqs = filter(lambda x: ".fq" in x, ls)

os.system("mkdir ../Split_reads")

for file in fqs:
        name = file.split('.')[0]
        os.system("paste - - - - - - - - < "+file+" |cut -f 1-4|tr '\t' '\n' >"+name+"_R1.fq")
        os.system("paste - - - - - - - - < "+file+" |cut -f 5-8|tr '\t' '\n' >"+name+"_R2.fq")

os.system("mv *R1.fq ../Split_reads/")
os.system("mv *R2.fq ../Split_reads/")
