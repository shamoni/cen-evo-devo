library(plyr)
library(ggplot2)

parseids <- function(x){
    X <-vector(mode="character",length=5)
    a <-strsplit(x,":")
    X[1] <-ifelse(substr(a[[1]][1],1,3) == "_R_",substr(x,4,nchar(a[[1]][1])),a[[1]][1])
    b <- strsplit(a[[1]][2],"-")
    X[2] <- b[[1]][1]
    X[3] <- b[[1]][2]
    X[4] <- b[[1]][3]
    X[5] <-ifelse(substr(a[[1]][1],1,3) == "_R_","2","1")
    return(X)
}

rep_pos_df <-function(x){
    df <-ldply(x,parseids)
    colnames(df) <-c("Contig","Start","Stop","Cluster","Strand")
    df$Start <- as.integer(df$Start)
    df$Stop <- as.integer(df$Stop)
    df <-df[with(df, order(Contig, Start)), ]
    d <-sapply(1:nrow(df), function(i) df$Start[i] - df$Stop[i-1])
    df$diff <- sapply(1:length(d), function(x) d[[x]][1])
    df$diff[!duplicated(df$Contig)] <- 0
    return(df)
}

parseclstrs <- function(df,feature){
    n <- nrow(df)
    ls <- character(0)
    for(i in 1:n){
        ls <- c(ls,df[i,feature])
    }
    return(ls)
}

parseclstrs_withgaps <- function(df,feature){
    n <- nrow(df)
    ls <- character(0)
    for(i in 1:n){
        if(df$diff[i] <= 175){ls <- c(ls,df[i,feature])}
        else{
            gap <- df$diff[i] %/% 175
            ls <- c(ls,df[i,feature],rep("7",gap))
        }
    }
    return(ls)
}

clstr_grid <- function(df,parse_fun,feature){
    out <- dlply(df,.(Contig),parse_fun,feature)
    dec <- names(sort(sapply(out,length)))
    out <- out[dec]
    ml <- max(sapply(out,length))
    out <-lapply(1:length(out),function(x) c(out[[x]],rep("0",ml-length(out[[x]]))))
    out.grid <- expand.grid(x=1:length(out[[1]]),y=1:length(out))
    out.grid$z <- unlist(out)
    return(out.grid)
}

plot_grid <- function(grid,feature){
    i<-ifelse(feature=="Cluster",1,2)
    col1 <- c("0"="#FFFFFF","1"="#FF4900FF","2"="#00A600FF","3"="#A0D600FF","4"="#FFDB00FF","5"="#00E5FFFF","6"="#0080FFFF","7"="#999999")
    col2 <-c("0"="#FFFFFF","1"="red","2"="blue","7"="#999999")
    col <-list(col1,col2)
    lab1 <- c("0"="","1"="Cluster1","2"="Cluster2","3"="Cluster3","4"="Cluster4","5"="Cluster5","6"="Cluster6","7"="Not CEN")
    lab2 <-c("0"="","1"="plus","2"="minus","7"="Not CEN180")
    lab <- list(lab1,lab2)
    p <- ggplot(grid,aes(x,rev(y),fill=z))
    p <- p + geom_raster() + scale_fill_manual(values=col[[i]],labels=lab[[i]])
    p <- p + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
    p <- p + theme(axis.ticks=element_blank(), axis.title=element_text(size=14))
    p <- p + labs(x="CEN180 repeats",y="")
    p <- p + theme(axis.text=element_blank())
    p <- p + guides(fill=guide_legend(title=""))
    return(p)
}

what_is_adjacent <- function(df){
    n <- nrow(df)
    ls <- character(0)
    for(i in 1:n){
        if(df$diff.Prev[i] <= 175 & !is.na(df$diff.Prev[i])){ls <- c(ls,df[i,"Prev"])}
        else{
            ls <- c(ls,"Not CEN180")
        }
        if(df$diff.Next[i] <= 175 & !is.na(df$diff.Next[i])){ls <- c(ls,df[i,"Next"])}
        else{
            ls <- c(ls,"Not CEN180")
        }
    }
    return(ls)
}

adj_counts_df <-function(x){
    df <-ldply(x,parseids)
    colnames(df) <-c("Contig","Start","Stop","Cluster","Strand")
    df$Start <- as.integer(df$Start)
    df$Stop <- as.integer(df$Stop)
    df <-df[with(df, order(Contig, Start)), ]
    d <-sapply(1:nrow(df), function(i) df$Start[i] - df$Stop[i-1])
    df$diff.Prev <- sapply(1:length(d), function(x) d[[x]][1])
    d <-sapply(1:nrow(df), function(i) df$Start[i+1] - df$Stop[i])
    df$diff.Next <- sapply(1:length(d), function(x) d[[x]][1])
    d <-sapply(1:nrow(df), function(i) df$Cluster[i-1])
    df$Prev <- sapply(1:length(d), function(x) d[[x]][1])
    d <-sapply(1:nrow(df), function(i) df$Cluster[i+1])
    df$Next <- sapply(1:length(d), function(x) d[[x]][1])
    i <- which(!duplicated(df$Contig))
    df$Prev[i] <- NA
    df$Next[i-1] <- NA
    out <- dlply(df,.(Cluster),what_is_adjacent)
    out.df <- lapply(as.character(1:6), function(x) as.data.frame(table(out[x]),stringsAsFactors=FALSE))
    freqs <-Reduce(function(...) merge(..., all=T,by="Var1"), out.df)
    freqs[is.na(freqs)] <- 0
    colnames(freqs) <-c("Adjacent","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6")
    return(freqs)
}

plot.adjacent <- function(df){
    custom <-c(heat.colors(10)[3],terrain.colors(10)[c(1,4)],heat.colors(10)[7],topo.colors(10)[c(4,3)],"#999999")
    p <-ggplot(df, aes(x=variable,y=value,fill=Adjacent)) + geom_bar(stat="identity")
    p <- p + scale_fill_manual(values=custom,labels=c(paste("Cluster",1:6,sep=""),"Not CEN180")) + theme_bw()
    p <- p + theme(axis.title=element_blank(),axis.text=element_text(size=12))
    p <- p + guides(fill=guide_legend(title="Adjacent Sequence"))
    return(p)
}

