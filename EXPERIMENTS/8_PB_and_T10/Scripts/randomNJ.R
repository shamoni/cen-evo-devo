library(phyclust)

seqdata <- read.fasta("q_CEN2_All.aln2")
findbest7_randomNJ <- find.best(seqdata$org,7,init.procedure=.init.procedure[-1],init.method=.init.method[3])
save(findbest7_randomNJ,file="findbest7_randomNJ.RData")
