library(parallel)

seqdata<- read.fasta("q_CEN2_All.aln2")

#Function to generate a list of 6 (500 random numbers) to sample the alignment
set.seed(123)
list.of.rnos <- function(alignment){
    rno_list <- list()
    for(i in 1:6){
        y <- sample(1:alignment$nseq,500, replace=FALSE)
        rno_list <- append(rno_list,list(y))
    }
    return(rno_list)
}

#Generate a list of nid matrixes based on rnos
sample <- list.of.rnos(seqdata)
nid_by_cluster <- list()
for(i in 1:length(sample)){
    temp <- seqdata$org[sample[[i]],]
    nid_by_cluster[[i]] <- temp
}

phyclust.result <-mclapply(nid_by_cluster, many_phyclust,ks=seq(2,10),mc.cores=6)

#Parses the output of phyclust.result into a list of phyclust objects and the dataframe
phyclust.clstrs <- lapply(1:length(phyclust.result), function(x) phyclust.result[[x]][[1]])
clstr.gram.input <- lapply(1:length(phyclust.result), function(x) phyclust.result[[x]][2])

#Generate a list of clustergrams
cg.list <- list()
for(i in 1:length(clstr.gram.input)){
    cg <- clustergram(clstr.gram.input[[i]][[1]])
    cg.list[i] <- list(cg)
}

save(phyclust.clstrs,clstr.gram.input,cg.list, file="2015-01-21-ClstrOptK2-10-PBandT10.RData")


