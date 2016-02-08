library(phyclust)

seqdata <- read.fasta("q_CEN2_All.aln2")
findbest7_randomMu <- find.best(seqdata$org,7,init.procedure=.init.procedure[-1],init.method=.init.method[1])
save(findbest7_randomMu,file="findbest7_randomMu.RData")
