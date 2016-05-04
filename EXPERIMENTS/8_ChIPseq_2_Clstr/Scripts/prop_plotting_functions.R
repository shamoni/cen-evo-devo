require(reshape)
require(ggplot2)
require(plyr)
require(RColorBrewer)

parse_clstr <- function(x){
    h <-strsplit(x,"-")[[1]][3]
}

parse_Chr <- function(x){
    h <-strsplit(x,":")[[1]][1]
    substr(h,nchar(h)-3,nchar(h))
}

cen.fraction <- function(df,meta){
    NChIP.summary <-adply(df, 1, transform, CEN = sum(Unique_CEN, Multi_CEN))
    metadata <- read.csv(meta, header=T)
    fileorder <- match(NChIP.summary$filename,metadata$LibraryID)
    ID.results <- metadata[fileorder,c(2,3)]
    mat <- data.matrix(NChIP.summary[,-c(1,3,4)])
    mat <- prop.table(mat,1)
    summary.prop <-cbind(Library=NChIP.summary$filename,ID.results,CEN=mat[,2])
    summary.means <- ddply(summary.prop,.(Sample,Antibody),summarise,mean=mean(CEN),sem=sd(CEN)/sqrt(length(CEN)))
    summary.means <- transform(summary.means, lower=mean-sem, upper=mean+sem)
    summary.means$Sample <- factor(summary.means$Sample,levels=c("WT (+/+)","LoCENH3 (-/-)","ZmCENH3 (-/-)"))
    ab <-c("INPUT","AtCENH3","LoCENH3","ZmCENH3")
    summary.means$Antibody <- factor(summary.means$Antibody,levels=ab)
    return(summary.means)
}

plot.cen.fraction <- function(df){
    p <-ggplot(data=df, aes(x=Sample,y=mean, fill=Sample)) + theme_bw()
    p <- p + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(),panel.grid=element_blank())
    p <- p + guides(fill=guide_legend(title="ChIP Sample"))
    p <- p + theme(strip.text.x = element_text(size = 8))
    p <- p + labs(y="CEN180 Reads/Total reads",x="")
    p <- p + theme(axis.title.y=element_text(size=10))
    p <- p + geom_bar(stat="identity") + geom_errorbar(aes(ymax=upper,ymin=lower), width=0.25)
    p <- p + scale_fill_manual(values=brewer.pal(9,"Blues")[c(5,7,9)])
    p <- p + facet_wrap(~ Antibody,nrow=1)
    return(p)
}

clstr.distribution <- function(df,meta,class){
    colnames(df) <- c("filename","Class","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6")
    metadata <- read.csv(meta, header=T)
    fileorder <- match(df$filename,metadata$LibraryID)
    ID.results <- metadata[fileorder,c(2,3)]
    mat <- data.matrix(df[,-c(1,2)])
    mat <- prop.table(mat,1)
    results.prop <-cbind(Library=df$filename,ID.results,Class=df$Class,as.data.frame(mat))
    results.melted <-melt(results.prop)
    to.omit <- c("SML_19_merged.fa","SML_38_merged.fa","SML_43_merged.fa")
    results.melted <-results.melted[-which(results.melted$Library %in% to.omit),]
    results <-results.melted[results.melted$Class == class,]
    results$group[which(results$Antibody == "INPUT")] <- "INPUT"
    results$group[which(results$Antibody == "AtCENH3")] <- "AtCENH3"
    results$group[which(results$Antibody == "LoCENH3")] <- "LoCENH3"
    results$group[which(results$Antibody == "ZmCENH3")] <- "ZmCENH3"
    results$tocolor <- results$group
    results$group <- factor(results$group,levels=c("INPUT","AtCENH3","LoCENH3","ZmCENH3"))
    results$tocolor <- factor(results$tocolor,levels=c("INPUT","AtCENH3","LoCENH3","ZmCENH3"))
    return(results)
}

plot.clstr.distribution <-function(df){
    custom <-c("black",brewer.pal(9,"Blues")[c(5,7,9)])
    p <- ggplot(data=df, aes(x=group,y=value)) + theme_bw()
    p <- p + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(),panel.grid=element_blank())
    p <- p + labs(y="Cluster-specific reads/All CEN180 reads")
    p <- p + theme(axis.title.y=element_text(size=10))
    p <- p + theme(legend.key=element_blank(), legend.title=element_blank())
    p <- p + geom_point(aes(fill=tocolor),colour="black",pch=21, size=5)
    p <- p + facet_wrap(~ variable,nrow=1) + scale_fill_manual(values=custom)
}




