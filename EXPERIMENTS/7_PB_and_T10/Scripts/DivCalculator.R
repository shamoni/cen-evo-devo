require(phyclust)

Hap.n <- function(haps){
    l <- length(haps)
    u <-length(sapply(1:l, function(x) haps[[x]][2]))
    return(u)   
}

compare.seq <- function(j,nid){
    seq.mat <- matrix(nid[j,],dim(nid)[1],dim(nid)[2],byrow=TRUE)
    diff.mat <- nid - seq.mat
    diff.mat[diff.mat !=0] <- TRUE
    diff.vect <- apply(diff.mat,1,sum)
    return(diff.vect)
}

Hap.list.by.cluster <- function(nid){
    n <-dim(nid)[1]
    skip <- numeric()
    haps <- list()
    for(i in 1:n){
        if(i %in% skip){
            next
        }
        else{
            differences <- compare.seq(i,nid)
            identical <- which(differences==0)
            skip <-c(skip,identical)
            row <- c(i,length(identical))
            haps <-c(haps,list(row))
        }
    }
    return(haps)
}

HaplotypeDiv <- function(Phyclst,X){
    n <- Phyclst$K
    Nid.by.n <- lapply(1:n, function(x) X[which(Phyclst$class.id ==x),])
    Out <-sapply(1:n, function(x) Hap.list.by.cluster(Nid.by.n[[x]]))
    Uniq <- sapply(1:n, function(x) Hap.n(Out[[x]]))
    return(Uniq)
}

Pi.calc <- function(nid){
    s <- sum(phyclust.edist(nid, edist.model = .edist.model[3]))
    n <-dim(nid)[1]
    l <-dim(nid)[2]
    pi <- 2*(1/n)*(1/(n-1))*s*(1/l)
    return(pi)
}

NucleotideDiv <- function(Phyclst,X){
    n <- Phyclst$K
    Nid.by.n <- lapply(1:n, function(x) X[which(Phyclst$class.id ==x),])
    Out <-sapply(1:n, function(x) Pi.calc(Nid.by.n[[x]]))
    return(Out)
}

