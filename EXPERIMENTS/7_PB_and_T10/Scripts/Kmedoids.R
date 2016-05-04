library(phyclust)

seqdata <- read.fasta("q_CEN2_All.aln2")
findbest7_Kmedoids <- find.best(seqdata$org,7,init.procedure=.init.procedure[-1],init.method=.init.method[5])
save(findbest7_Kmedoids,file="findbest7_Kmedoids.RData")
