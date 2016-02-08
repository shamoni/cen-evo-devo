#! /usr/bin/env python

import os, sys, math

#Generate a list with .fq files in current folder
ls = os.listdir(os.getcwd())
fqs = filter(lambda x: ".fq" in x, ls)

os.system("mkdir ../After_CDHIT-DUP")

for file in fqs:
        name = "_".join(file.split('_')[:-1])
        os.system("~/bin/cdhit-master/cd-hit-auxtools/cd-hit-dup -i "+file+" -o "+name+"_cdhitdup -m false  -u 0")

os.system("rm *2.clstr")
os.system("mv *cdhitdup ../After_CDHIT-DUP/")
os.system("mv *.clstr ../After_CDHIT-DUP/")

