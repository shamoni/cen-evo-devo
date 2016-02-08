require(reshape)
require(ggplot2)
require(plyr)

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
    to.omit <- c("SML_12_merged.fa")
    summary.prop <-summary.prop[-which(summary.prop$Library %in% to.omit),]
    summary.means <- ddply(summary.prop,.(Sample,Antibody),summarise,mean=mean(CEN),sem=sd(CEN)/sqrt(length(CEN)))
    summary.means <- transform(summary.means, lower=mean-sem, upper=mean+sem)
    summary.means$Sample <- factor(summary.means$Sample,levels=c("WT (+/+)","LoCENH3 (-/-)"))
    ab <-c("INPUT","H3K36me2","H3K27me3","H3K9me2","H3K27me1","AtCENH3","LoCENH3")
    summary.means$Antibody <- factor(summary.means$Antibody,levels=ab)
    return(summary.means)
}

plot.cen.fraction <- function(df){
    p <-ggplot(data=df, aes(x=Sample,y=mean, fill=Sample)) + theme_bw()
    p <- p + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(),panel.grid=element_blank())
    p <- p + theme(legend.title=element_blank())
    p <- p + theme(strip.text.x = element_text(size = 8))
    p <- p + geom_bar(stat="identity") + geom_errorbar(aes(ymax=upper,ymin=lower), width=0.25)
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
    to.omit <- c("SML_12_merged.fa")
    results.melted <-results.melted[-which(results.melted$Library %in% to.omit),]
    results.unique <-results.melted[results.melted$Class == class,]
    results.unique$group <- "Other marks"
    results.unique$group[which(results.unique$Antibody == "INPUT")] <- "INPUT"
    results.unique$group[which(results.unique$Antibody == "AtCENH3")] <- "AtCENH3"
    results.unique$group[which(results.unique$Antibody == "LoCENH3")] <- "LoCENH3"
    results.unique$tocolor <- results.unique$group
    results.unique$tocolor[which(results.unique$Library %in% c("SML_19_merged.fa","SML_26_merged.fa"))]<-"Control"
    results.unique$group <- factor(results.unique$group,levels=c("INPUT","Other marks","AtCENH3","LoCENH3"))
    results.unique$tocolor <- factor(results.unique$tocolor,levels=c("INPUT","Other marks","AtCENH3","LoCENH3","Control"))
    return(results.unique)
}

plot.clstr.distribution <-function(df){
    custom <-c(gray(1:10/10)[c(1,5)],heat.colors(10)[c(1,5)],cm.colors(10)[1])
    p <- ggplot(data=df, aes(x=group,y=value)) + theme_bw()
    p <- p + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(),panel.grid=element_blank())
    p <- p + theme(legend.key=element_blank(), legend.title=element_blank())
    p <- p + geom_point(aes(fill=tocolor),colour="black",pch=21, size=5)
    p <- p + facet_wrap(~ variable,nrow=1) + scale_fill_manual(values=custom)
}




