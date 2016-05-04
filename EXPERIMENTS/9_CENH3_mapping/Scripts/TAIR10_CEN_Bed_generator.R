source("../../8_ChIPseq_2_Clstr/Scripts/clstr_pos_functions.R")

library(phyclust)

seqdata <- read.fasta("../../8_ChIPseq_2_Clstr/q_CEN2_All_K6.aln2")
ids <- seqdata$seqname
ids.T10 <- ids[grepl("Chr",ids)]
T10.df <-rep_pos_df(ids.T10)

Bed.df <- T10.df[,c(1,2,3)]
names(Bed.df) <- c("chrom","chromStart","chromEnd")
Bed.df$name <- paste("Cluster",T10.df$Cluster,sep="")
Bed.df$score <- 0
Bed.df$strand <- ifelse(T10.df$Strand == 1,"+","-")
Bed.df$thickStart <- Bed.df$chromStart
Bed.df$thickEnd <- Bed.df$chromEnd
Bed.df$itemRgb <- T10.df$Cluster

Bed.df$itemRgb[Bed.df$itemRgb==1] <- "239,77,34"
Bed.df$itemRgb[Bed.df$itemRgb==2] <- "15,165,74"
Bed.df$itemRgb[Bed.df$itemRgb==3] <- "160,204,58"
Bed.df$itemRgb[Bed.df$itemRgb==4] <- "253,218,0"
Bed.df$itemRgb[Bed.df$itemRgb==5] <- "94,202,230"
Bed.df$itemRgb[Bed.df$itemRgb==6] <- "67,119,188"

write.table(Bed.df, file="../Cluster.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)