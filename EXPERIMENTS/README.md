This directory contains folders with either defunct or ongoing experiments, for e.g. trying out new analysis methods, parameter optimization, developing visualization protocols. The documentation for these methods is more of an overview rather than a step-by-step protocol. Fruitful finalized strategies are migrated to the RESULTS directory.

## Estimate CEN fraction by mapping reads to an arbitrary reference using Bowtie
	/EXPERIMENTS/1_Bowtie_mapto_180X2

I generated an arbitrary reference genome from a tandem duplicate of the centromere repeat sequence i.e. 180X2. I then built the bowtie2 index for this genome. I then used the script bowtie-very-sensitive.py to process all ChIPseq reads in a folder through the following pipeline: map to 180X2 reference using the --very-sensitive parameters, use smalls to generate the bam and bay files and capture the standard out. The standard out from each mapping is moved to a separate folder called STDOUT and manually processed to generate a summary data table.

## Test the Henikoff strategy for identifying functional centromere sequences by CDHIT clustering
	/EXPERIMENTS/2_SH_CDHIT_Clstrs   
	
The strategy is implemented in [Henikoff2015]. In this publication the default cd-hit clustering threshold of 0.9 was used, I first wanted to test if there was an optimal threshold for clustering the ChIPseq reads. I ran a script cdhit-est-g1r1.py to cluster SML_17 ChIPseq reads at thresholds ranging from 1.0 to 0.9 at -0.01 increments. I also implemented additional parameters: g 1 (not greedy instead find the best cluster) and r 1 to search the reverse strand. I used reads from After_CDHIT-DUP as the starting point. This analysis did not reveal an inflection point suggesting that a threshold of 0.9 is suitable.

I then used cd-hit-est to cluster the following libraries: SML_11_dup.fa; SML_17_dup.fa; SML_18_dup.fa; SML_19_dup.fa

	module load cdhit
	cd-hit-est -i SML_11_dup.fa -o SML_11_.90 -n 10 -c 0.9 -d 0 -M 0 -T 16

The 0.9 cluster references from SML_17 (AtCENH3 ChIP) were used to generate a BWA indexed genome in the folder:17_0.90_cluster_ref
   
   	module load bwa
	bwa index -p 17_.90 -a is SML_17_.90.fa

The goal is to identify the clusters representing the most abundant reads in a ChIPseq library. To do so simply map reads to these clusters as references using bra -n 1 parameter (i.e. keeping only one hit)

	bwa aln -t 10 17_.90 SML_17_merged.fq > SML_17.sai
	bwa samse -n 10 17_.90 SML_17.sai SML_17_merged.fq > SML_17.sam
	samtools view -b -S -F 4 SML_17.sam > SML_17.bam
	samtools sort SML_17.bam SML_17_sorted
	samtools index SML_17_sorted.bam
	samtools idxstats SML_17_sorted.bam > 17idx.out
	sort -nrk3 17idx.out > 17idx_sorted.out

Used the python script rank-extract.py to grab ids of clusters with no. of reads mapped is > threshold. Then using the list of ids generate a fasta file with sequences of the select top ranking clusters.

## Identify PacBio contigs with centromere repeats using BLAST
	/EXPERIMENTS/3_PacBio_BLAST
The 13 SML_17 clusters which had > 50,000 SML_17 ChIPseq reads were selected as the reference sequences. A BLAST database was generated using these sequences and all downloaded PacBio contigs were used as the query.

	makeblastdb -in 17_top50K_clstr_ref.fa -dbtype nucl -out 17_top50K
	blastn -task megablast -db 17_top50K -query JSAD01.1.fa -out q_JSAD01_db_17top50K_blast -perc_identity 90 -word_size 50 -outfmt “6 qseqid qlen sseqid qstart qend sstrand length pident mismatch gapopen evalue bitscore”
	cut -f1 q_JSAD01_db_17top50K_blast | uniq > JSAD01_17top50K_ids
	sort -k1,1r -k4,4n -k8 q_JSAD01_db_17top50K_blast > q_JSAD01_db_17top50K_blast_sorted

## Analysis of PacBio contig006
	/EXPERIMENTS/4_Contig006_PacBio
This was the top PacBio contig that the previous BLAST analysis spit out. I made a folder with contig006 as the BWA indexed reference genome and mapped reads from ChIPseq libraries SML_11, SML_14, SML_17, SML_18 and SML_19 to it.

	bwa aln -t 5 contig006/contig6 /Merged_reads/SML_14_merged.fq > 14_mapto_contig6.sai
	bwa samse -n 10 contig006/contig6 11_mapto_contig6.sai /Merged_reads/SML_14_merged.fq > 14_mapto_contig6.sam
	module load samtools/default
	samtools view -b -S -F 4 14_mapto_contig6.sam > 14_mapto_contig6.bam
	samtools sort 14_mapto_contig6.bam 14_mapto_contig6_sorted
	rm 14_mapto_contig6.bam
	samtools index 14_mapto_contig6_sorted.bam
	~/bin/IGVTools/igvtools count --includeDuplicates -w 1 14_mapto_contig6_sorted.bam 14_mapto_contig6.tdf contig006/contig006.fa
	
The mapping showed a really pronounced phasing pattern at the 5' end of contig006. This was also the region where the majority of the reads from all the libraries were mapping. To test if this was the effect of its position along the contig I repeated the mapping experiment but this time I using a contig006 reference in which this "regionA" was masked.

	module load bedtools2/2.21.0
	bedtools maskfasta -fi ../contig006.fa -bed regionA.bed -fo contig6_regionA_masked.fa
	module load bwa
	bwa index -p contig6_regionA_masked -a is contig6_regionA_masked.fa
	....

This mapping to the rest of the contain remained unchanged after masking the 5 'end, suggesting that it is not an effect of its position being first on the contain. I then wanted to discover all repeats in contig006 and compare them to the sequences in regionA. To do so I ran contig006 through [Tandem Repeat Finder](https://tandem.bu.edu/trf/trf.basic.submit.html) which identified five *indices* with period sizes ranging from 176-178 bp. I grabbed the consensus sequences from each of these indices and made a bed file of the five indices to indicate regions contig006 with tandem CEN repeats. 

Contig | Start | Stop | Region
--- | --- | --- | --- |
JSAD01000006.1|	0 | 7915 | A |
JSAD01000006.1|	9496	 | 11094 | B 
JSAD01000006.1|	27096 | 35798 |	 C 
JSAD01000006.1|	58934 | 78638 |	 D 
JSAD01000006.1|	94407 | 100136 | E 


I made tandem duplicates from each of the consensus sequences i.e. each sequence represented two consecutive repeats. I then aligned them to each other along with sequence from one of the phased CEN allele rotected by CENH3 ChIPseq (Block A mapping). This file is called contig6_TRF_consensus_inregister.fa. These in register consensus repeats were made into a blast database and queried with contig006. The output file was sorted and locations of highest-scoring non-overlapping repeats was extracted in bed format using the script extract-repeats.py. A total of 233 repeats were identified in contig006 which has a length of 100137 bp (~ 41% coverage). The number of non-overlapping repeats identified in contig006 by the top50K clusters as the BLAST reference was only 34. This is likely because the top50K cluster references from AtCENH3 ChIP are highly homogeneous and only identify repeats in regionA.

	module load bedtools2/2.21.0
	bedtools getfasta -name -fi contig006.fa -bed contig006_TRFinregister_repeats.bed -fo contig006_repeats.fa

## Identify and align centromere repeats in PacBio contigs using LASTZ
	/EXPERIMENTS/5_PacBio_LASTZ

I copied the contig6_TRF_consensus_inregister.fa file from 4_Contig006_PacBio/contig006_BLAST_db/ folder and renamed the IDs CEN1, CEN2...CEN6. I then queried all PacBio contigs using these six repeats as query sequences. I parsed the output using the parse_lastz.py

	lastz ../../DATA/2015-PacBio/JSAD01.1.fa[nameparse=darkspace][multiple] contig6_TRF_consensus_inregister.fa --coverage=90 --format=general:score,name1,strand1,size1,start1,end1,name2,strand2,identity,length1,align1 > q_contig6TRF_PacBio.csv
	../../PROTOCOLS/parse_lastz.py q_contig6TRF_PacBio.csv -o -s -p
	Recovered 28013 records

After checking the output from all six query sequences I decided to move forward with only CEN2 as query.
	
	lastz ~/DATA/2015-PacBio/JSAD01.1.fa[nameparse=darkspace][multiple] CEN2.fa --coverage=90 --format=general:score,name1,strand1,size1,start1,end1,name2,strand2,identity,length1,align1 > q_CEN2_PacBio.csv
	~/PROTOCOLS/parse_lastz.py -r -g 14 q_CEN2_PacBio.csv 
	Recovered 4780 records
	Recovered 4725 records with less than 14 gaps in target sequence
	Recovered 4547 records of lengths 165-180 bp
	
I used clustalo to align the repeats and msa_refiner.py to refine the alignment.

	~/PROTOCOLS/msa_refiner.py q_CEN2_PacBio.aln
	Recovered msa of lengh 240 with 4547 aligned sequences
	Kept 4493 aligned sequences after deleting sequences with insertions shared by < 20 other sequences
	Trimmed msa has length 202 with 4493 aligned sequences

## Use LASTZ to extract repeats from the Arabidopsis Genome Assembly i.e. TAIR10
	/EXPERIMENTS/6_TAIR10_LASTZ
	
	lastz /isner/share/genomes/thaliana/TAIR10/Genomic/TAIR10_all.fas[nameparse=darkspace][multiple] CEN2.fa --coverage=90 --ambiguous=iupac --format=general:score,name1,strand1,size1,start1,end1,name2,strand2,identity,length1,align1  > q_CEN2_TAIR10.csv
	../../PROTOCOLS/parse_lastz.py q_CEN2_TAIR10.csv -o -s -p
	Recovered 2550 records
	
	../../PROTOCOLS/parse_lastz.py -r -g 14 q_CEN2_TAIR10.csv
	Recovered 2550 records
	Recovered 2541 records with less than 14 gaps in target sequence
	Recovered 2458 records of lengths 165-180 bp
	clustalo -i q_CEN2_TAIR10.fa -o q_CEN2_TAIR10.aln
	~/PROTOCOLS/msa_refiner.py q_CEN2_TAIR10.aln
	Recovered msa of lengh 203 with 2458 aligned sequences
	Kept 2399 aligned sequences after deleting sequences with insertions shared by < 20 other sequences
	Trimmed msa has length 189 with 2399 aligned sequences
	
## Generate clusters by analyzing all repeats identified i.e. both PacBio and TAIR10
	/EXPERIMENTS/7_PB_and_T10
	
	cp ../5_PacBio_LASTZ/q_CEN2_PacBio.fa .
	cp ../7_TAIR10_LASTZ/q_CEN2_TAIR10.fa .
	cat q_CEN2_PacBio.fa q_CEN2_TAIR10.fa > q_CEN2_All.fa
	rm q_CEN2_PacBio.fa
	rm q_CEN2_TAIR10.fa
	clustalo -i q_CEN2_All.fa -o q_CEN2_All.aln
	../../PROTOCOLS/msa_refiner.py -t 50 q_CEN2_All.aln 
	Recovered msa of lengh 234 with 7005 aligned sequences
	Kept 6868 aligned sequences after deleting sequences with insertions shared by < 50 other sequences
	Trimmed msa has length 204 with 6868 aligned sequences
	
To optimize the clustering I attempted cluster a subsets of the data from K=2 to K=10 using the find.best() phyclust function.

	R
	source("../../PROTOCOLS/my_clustergram.R")
	source("Clgram_K2K10_run.R")

To generate the best possible clusters I tested phyclust find.best algorithm with different init.procedure values and three init.method strategies: K-medoids, randomNJ, randomMu for both K=6 and K=7

	R
	source("Kmedoids.R")
	source("randomNJ.R")
	source("randomMu.R")

To write a fast file in which the cluster id is appended to the name of the sequence

	R
	library(phyclust)
	seqdata <- read.fasta("q_CEN2_All.aln2")
	load("findbest6_Kmedoids.RData")
	write.fasta(seqdata=seqdata$org,"q_CEN2_All_K6.aln2",classid=findbest6_Kmedoids$class.id,seqname= seqdata$seqname)
	mv q_CEN2_All_K6.aln2 ../9_ChIPseq_2_Clstr

## Assign ChIPseq reads to CEN clusters and analyze distribution of reads across clusters
	/EXPERIMENTS/8_ChIPseq_2_Clstr
	~/PROTOCOLS/parse_in2_clstrs.py q_CEN2_All_K6.aln2

	~/cen-evo-devo/PROTOCOLS/Seq_2_Clstr.py -Fc q_CEN2_All_K6_Clusters/ -Fr ~/cen-evo-devo/DATA/Select_Data/ -O NChIP-K30 -K 30

	Processing Kmers from 6 Clusters
	Identifying Kmers of size 30 shared by at least 100 sequences within the cluster

	q_CEN2_All_K6.aln2_Cluster1.fa; No. of Kmers 138; No. of Signature Kmers 7
	q_CEN2_All_K6.aln2_Cluster2.fa; No. of Kmers 206; No. of Signature Kmers 166
	q_CEN2_All_K6.aln2_Cluster3.fa; No. of Kmers 387; No. of Signature Kmers 304
	q_CEN2_All_K6.aln2_Cluster4.fa; No. of Kmers 153; No. of Signature Kmers 21
	q_CEN2_All_K6.aln2_Cluster5.fa; No. of Kmers 99; No. of Signature Kmers 99
	q_CEN2_All_K6.aln2_Cluster6.fa; No. of Kmers 350; No. of Signature Kmers 350

	~/cen-evo-devo/PROTOCOLS/Seq_2_Clstr.py -Fc q_CEN2_All_K6_Clusters/ -Fr ~/cen-evo-devo/DATA/Select_Data/ -O NChIP-K25 -K 25

	Processing Kmers from 6 Clusters
	Identifying Kmers of size 25 shared by at least 100 sequences within the cluster

	q_CEN2_All_K6.aln2_Cluster1.fa; No. of Kmers 161; No. of Signature Kmers 18
	q_CEN2_All_K6.aln2_Cluster2.fa; No. of Kmers 249; No. of Signature Kmers 182
	q_CEN2_All_K6.aln2_Cluster3.fa; No. of Kmers 447; No. of Signature Kmers 335
	q_CEN2_All_K6.aln2_Cluster4.fa; No. of Kmers 161; No. of Signature Kmers 16
	q_CEN2_All_K6.aln2_Cluster5.fa; No. of Kmers 127; No. of Signature Kmers 127
	q_CEN2_All_K6.aln2_Cluster6.fa; No. of Kmers 434; No. of Signature Kmers 434

	~/cen-evo-devo/PROTOCOLS/Seq_2_Clstr.py -Fc q_CEN2_All_K6_Clusters/ -Fr ~/cen-evo-devo/DATA/Select_Data/ -O NChIP-K20 -K 20

	Processing Kmers from 6 Clusters
	Identifying Kmers of size 20 shared by at least 100 sequences within the cluster

	q_CEN2_All_K6.aln2_Cluster1.fa; No. of Kmers 176; No. of Signature Kmers 22
	q_CEN2_All_K6.aln2_Cluster2.fa; No. of Kmers 296; No. of Signature Kmers 182
	q_CEN2_All_K6.aln2_Cluster3.fa; No. of Kmers 477; No. of Signature Kmers 329
	q_CEN2_All_K6.aln2_Cluster4.fa; No. of Kmers 169; No. of Signature Kmers 14
	q_CEN2_All_K6.aln2_Cluster5.fa; No. of Kmers 144; No. of Signature Kmers 138
	q_CEN2_All_K6.aln2_Cluster6.fa; No. of Kmers 480; No. of Signature Kmers 475

## Map ChIPseq reads to the TAIR10 reference genome and analyze peaks and correlation between binding profiles
	/EXPERIMENTS/9_CENH3_mapping
	.../../PROTOCOLS/bwa-mapping.py -f ../../DATA/Select_Data/ -d T10_index/T10 -o BWA_mapping
	../../PROTOCOLS/macs_peakfinding.sh BWA_mapping/SML_17_mapto_T10.bam BWA_mapping/SML_19_mapto_T10.bam AtCENH3_rep1
	../../PROTOCOLS/macs_peakfinding.sh BWA_mapping/SML_18_mapto_T10.bam BWA_mapping/SML_19_mapto_T10.bam AtCENH3_rep2
	../../PROTOCOLS/macs_peakfinding.sh BWA_mapping/SML_23_mapto_T10.bam BWA_mapping/SML_43_mapto_T10.bam LoCENH3_rep1
	../../PROTOCOLS/macs_peakfinding.sh BWA_mapping/SML_25_mapto_T10.bam BWA_mapping/SML_43_mapto_T10.bam LoCENH3_rep2
	../../PROTOCOLS/macs_peakfinding.sh BWA_mapping/SML_39_mapto_T10.bam BWA_mapping/SML_38_mapto_T10.bam ZmCENH3_rep1
	../../PROTOCOLS/macs_peakfinding.sh BWA_mapping/SML_40_mapto_T10.bam BWA_mapping/SML_38_mapto_T10.bam LoCENH3_rep2
	
	mkdir DiffBind_Analysis
Generate the file cenh3_metadata.csv
	
	R
	library("DiffBind")	
Construct a new dba object from sample sheet.

	cenh3 <-dba(sampleSheet="cenh3_metadata.csv")	
Counts reads in binding site intervals, by default only includes peaks that are present in at least two peaksets. The binding affinity matrix represents scores calculated by TMM normalization (using edgeR), using ChIP read counts minus Control read counts and full library size.

	cenh3_affinity <- dba.count(cenh3)
	save(cenh3_affinity, file="cenh3_affinity.RData")
	
	library(corrplot)
	corvals <- dba.plotHeatmap(cenh3_affinity)
	pdf("corrplot.pdf")
	corrplot(corvals,method="shade",shade.col=NA,tl.col="black",tl.srt=45,addCoef.col="black",order="AOE")
	dev.off()
	
Generate a bed file for the CEN180 repeats colored by cluster, this will be useful for comparing along with the CENH3 Fold Enrichment bigwig (bw) files on IGV

	cd Scripts/
	R
	source("TAIR10_CEN_Bed_generator.R")

	
[Henikoff2015]: http://www.ncbi.nlm.nih.gov/pubmed/25927077


