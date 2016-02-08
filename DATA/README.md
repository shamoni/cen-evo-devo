This directory contains all the data, both downloaded and self-generated, related to the *A. thaliana* centromere project. The self-generated data are native histone ChIPseq experiments. Subdirectories within each ChIPseq directory contains reads in different stages of the processing pipeline.

METADATA
----------------
The following files/links correspond to the metadata for each of the datasets

* 2015-03-NChIP/2015-03-NChIP-metadata.csv
* 2015-08-NChIP/2015-08-NChIP-metadata.csv
* 2015-BAC   
[Kumekawa2000] and [Kumekawa2001] did fine mapping of BACs to the centromeres of *A. thaliana* chromosome 5 and 4. I do not remember from exactly where I downloaded the BACs but here are links to the most likely places:   
	<ftp://ftp.arabidopsis.org/home/tair/Sequences/clones/>
	<ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/>
	<ftp://ftp.ncbi.nih.gov/repository/clone/reports/Arabidopsis_thaliana/>
* 2015-PacBio<br/>A fast file of PacBio contigs from [Kim2014] was downloaded from this link:   
	<http://www.ncbi.nlm.nih.gov/Traces/wgs/?val=JSAD01&page=1&display=contigs&search=JSAD01000110>

ChIPseq PROCESSING PIPELINE
---------------------------
The scripts used for batch processing of the reads is in the folder Data/Processing_Scripts
	#Use allprep for quality filter and demultiplex the lane
	/isner/share/scripts/allprep-10-current.py -b Barcode_ChIP_Aug2015.txt -f SM-ChIP-Aug2015_S74_L007_R1_001.fastq -r SM-ChIP-Aug2015_S74_L007_R3_001.fastq -i SM-ChIP-Aug2015_S74_L007_R2_001.fastq -m

	#Move the raw reads into their own folder, including the barcode data file and gzip to save space.
	mkdir Raw_reads/ 
	mv *.fastq Raw_reads/ 
	gzip Raw_reads/*.fastq 
	mv Barcode_ChIP_Aug2015.txt Raw_reads/

	#Move the filtered reads into their own folder and split the files into forward and reverse reads. This format is needed for the SeqPrep. 
	mkdir Reads_after_QC 
	mv *.fq Reads_after_QC/ 
	cd Reads_after_QC/ 
	split_interleaved_fq.py

	#Generate merged reads, tar compress the split reads, unzip merged fq reads, generate corresponding .fa files for merged reads.
	cd ../Split_reads
	merge_reads.py
	tar -czvf split_reads.tar.gz *.fq
	rm *.fq
	cd ../Merged_reads
	gunzip *.gz
	fq-2-fa.py
	module load fastqc #Module java-jdk1.8.0_05-static loaded.#Module fastqc-v0.11.2-static loaded.
	mkdir fastQC
	fastqc -o fastQC/ *.fq

	#Estimate level of sequence duplicates using CDHIT-DUP and generate reads with unique sequences in folder After_CDHITDUP
	cdhitdup_fq.py
	cd After_CDHIT-DUP/
	for file in *dup; do mv $file $file.fq; done
	fq-2-fa.py

[Kumekawa2000]: http://www.ncbi.nlm.nih.gov/pubmed/11214966
[Kumekawa2001]: http://www.ncbi.nlm.nih.gov/pubmed/11853315
[Kim2014]: http://www.nature.com/articles/sdata201445
