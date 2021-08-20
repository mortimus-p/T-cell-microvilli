
############Data preprocessing############
#Read  matrix files (after transcript mapping and counting)
get_files<-list.files()
get_sample<-read.delim(file=get_files[1],header = TRUE, sep="\t",dec = ".")


dataset<-as.data.frame(get_sample)
rownames(dataset)<-dataset$Identifier
dataset[,1]<-dataset[,2]

for ( i in 2:length(get_files)){
  get_sample<-read.delim(file=get_files[i],header = TRUE, sep="\t",dec = ".")
  
  
  dataset[,i]<- get_sample[,2]
  
}


colnames(dataset)<-sub("\\.txt", "", get_files)


# change transcript annotation to gene symbols
library(AnnotationHub)
library(fgsea)
library("GSEABase")
library(qusage)

ah = AnnotationHub()
ah <- subset(ah, species == "Homo sapiens")
query(ah, "GRCh38")



edb <- ah[["AH75011"]] #  database for human
columns(edb)
geneIDs <- rownames(dataset)
get_geneID<-select(edb, keys=geneIDs, columns= c("GENEID","DESCRIPTION", "GENENAME"),keytype="GENEID")
toreplace<-geneIDs
toreplace[which(toreplace %in% get_geneID$GENEID)]<-get_geneID$GENENAME
rownames(dataset)<-make.unique(toreplace)


# read sample info from file. This file contains sample names and conditions (surface and stimulation)
get_info<-read.delim("sample_info.txt", header=TRUE, sep="\t")
get_info$Group<-factor(paste(get_info$stim,get_info$surface, sep="."))
design<-model.matrix(~0+get_info$Group)
colnames(design) <- levels(get_info$Group)

y <- DGEList(counts=dataset, group=c(1,1,1,2,2,2,3,3,3,4,4,4))
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

#TMM normalization (to account for different sequencing depths)
y <- calcNormFactors(y)
y$samples


y <- estimateDisp(y, design, robust=TRUE) 
y$common.dispersion



fit <- glmQLFit(y, design, robust=TRUE) 



pdf(file=paste("MDS.pdf",sep=""),height = 10, width=10)
plotMDS(y, col=rep(1:4, each=3))
dev.off()


save(y, file="y.Robj")
save(anov, file="anov.Robj")
save(design, file="design.Robj")
save(fit, file="fit.Robj")



con<-makeContrasts(
  porous_US_flat_Ab = US.porous - AB.flat,
  porous_flat_US = US.porous - US.flat,
  AB_porous_US_flat = AB.porous - US.flat,
  porous_flat_AB = AB.porous - AB.flat,
  AB_US_porous = AB.porous - US.porous,
  AB_US_flat = AB.flat - US.flat,
  porous_flat = (AB.porous + US.porous) - (AB.flat + US.flat),
  AB_US = (AB.flat+ AB.porous) - (US.flat + US.porous),
  AB_only_porous = (AB.porous - US.porous) - (AB.flat -US.flat),
  porousAb_only = AB.porous - (AB.flat+ US.porous+US.flat)/3,
  levels=design
)



res_porous_US_flat_Ab<-glmQLFTest(fit, contrast=con[,"porous_US_flat_Ab"])
res_porous_flat_US<-glmQLFTest(fit, contrast=con[,"porous_flat_US"])
res_porous_flat_AB<-glmQLFTest(fit, contrast=con[,"porous_flat_AB"])
res_AB_porous_US_flat<-glmQLFTest(fit, contrast=con[,"AB_porous_US_flat"])
res_AB_US_porous<-glmQLFTest(fit, contrast=con[,"AB_US_porous"])
res_AB_US_flat<-glmQLFTest(fit, contrast=con[,"AB_US_flat"])
res_AB_US_flat<-glmQLFTest(fit, contrast=con[,"AB_US_flat"])
res_porous<-glmQLFTest(fit, contrast=con[,"porous_flat"])
res_AB_US<-glmQLFTest(fit, contrast=con[,"AB_US"])
res_AB_only_porous<-glmQLFTest(fit, contrast=con[,"AB_only_porous"])
res_porousAb_only<-glmQLFTest(fit, contrast=con[,"porousAb_only"])


logCPM <- cpm(y, prior.count=10, log=TRUE)




is.de<-decideTestsDGE(res_porous_US_flat_Ab, lfc=0.5)
summary(is.de)
plotMD(res_porous_US_flat_Ab, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_porous_US_flat_Ab$table
DE_genes$symbol<-res_porous_US_flat_Ab$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]
write.csv(DE_genes,
          file="DE_porous_US_flat_Ab.csv",
          quote=FALSE)

is.de<-decideTestsDGE(res_porous_flat_US, lfc=0.5)
summary(is.de)
plotMD(res_porous_flat_US, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_porous_flat_US$table
DE_genes$symbol<-res_porous_flat_US$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]
write.csv(DE_genes,
          file="DE_genes_porous_flat_US.csv",
          quote=FALSE)


is.de<-decideTestsDGE(res_porous_flat_AB, lfc=0.5)
summary(is.de)
plotMD(res_porous_flat_AB, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_porous_flat_AB$table
DE_genes$symbol<-res_porous_flat_AB$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]
write.csv(DE_genes,
          file="DE_genes_porous_flat_AB.csv",
          quote=FALSE)


is.de<-decideTestsDGE(res_AB_US_porous, lfc=0.5)
summary(is.de)
plotMD(res_AB_US_porous, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_AB_US_porous$table
DE_genes$symbol<-res_AB_US_porous$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]
write.csv(DE_genes,
          file="DE_genes_AB_US_porous.csv",
          quote=FALSE)




is.de<-decideTestsDGE(res_AB_US_flat, lfc=0.5)
summary(is.de)
plotMD(res_AB_US_flat, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_AB_US_flat$table
DE_genes$symbol<-res_AB_US_flat$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]
write.csv(DE_genes,
          file="DE_genes_AB_US_flat.csv",
          quote=FALSE)



is.de<-decideTestsDGE(res_AB_porous_US_flat, lfc=0.5)
summary(is.de)
plotMD(res_AB_porous_US_flat, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_AB_porous_US_flat$table
DE_genes$symbol<-res_AB_porous_US_flat$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]
write.csv(DE_genes,
          file="DE_genes_AB_porous_US_flat.csv",
          quote=FALSE)


is.de<-decideTestsDGE(res_porous, lfc=0.5)
summary(is.de)
plotMD(res_porous, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_porous$table
DE_genes$symbol<-res_porous$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]
write.csv(DE_genes,
          file="DE_genes_porous_flat.csv",
          quote=FALSE)



is.de<-decideTestsDGE(res_AB_US, lfc=0.5)
summary(is.de)
plotMD(res_AB_US, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_AB_US$table
DE_genes$symbol<-res_AB_US$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]

write.csv(DE_genes,
          file="DE_genes_AB_US.csv",
          quote=FALSE)


is.de<-decideTestsDGE(res_AB_only_porous, lfc=0.5)
summary(is.de)
plotMD(res_AB_only_porous, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_AB_only_porous$table
DE_genes$symbol<-res_AB_only_porous$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]

write.csv(DE_genes,
          file="DE_genes_AB_only_porous.csv",
          quote=FALSE)





is.de<-decideTestsDGE(res_porousAb_only, lfc=0.5)
summary(is.de)
plotMD(res_porousAb_only, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
DE_genes<-res_porousAb_only$table
DE_genes$symbol<-res_porousAb_only$genes
DE_genes<-DE_genes[which(is.de@.Data!=0),]

write.csv(DE_genes,
          file="DE_genes_porousAb_only.csv",
          quote=FALSE)




############Perform GO term enrichment on the different lists############

library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)
library(Rgraphviz)
library("GSEABase")

#read files containing differentially expressed genes between different groups
getfiles<-list.files("DE_genes")
for (i in 1:length(getfiles)){
  DE_genes<-read.csv(paste("DE_genes/",getfiles[i],sep=""))  
  
  
  filt_up<-subset(DE_genes,DE_genes$logFC>0)
  filt_down<-subset(DE_genes,DE_genes$logFC<0)
  
  
  gene.df <- bitr(filt_up$X, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  
  ego <- enrichGO(gene         = gene.df$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  # keytype       = 'SYMBOL',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable =TRUE)

  ego2 <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   # keytype       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable = TRUE)
  
  
  ego3 <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   # keytype       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable=TRUE)
  
  
  gene.df <- bitr(filt_down$X, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  
  ego4 <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable =TRUE)
  
  
  ego5 <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable = TRUE)
  
  
  ego6 <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable=TRUE)
  
  
  pdf(paste("Path_en_",getfiles[i],"_up_down.pdf",sep=""))
  print(dotplot(ego,showCategory=20))
  print(dotplot(ego2,showCategory=20))
  print(dotplot(ego3,showCategory=20))
  
  print(dotplot(ego4,showCategory=20))
  print(dotplot(ego5,showCategory=20))
  print(dotplot(ego6,showCategory=20))
  dev.off()
  
  
  write.csv(ego, file=paste("CC_up",getfiles[i], sep=""))
  write.csv(ego2, file=paste("BP_up",getfiles[i], sep=""))
  write.csv(ego3, file=paste("MF_up",getfiles[i], sep=""))
  write.csv(ego4, file=paste("CC_down",getfiles[i],  sep=""))
  write.csv(ego5, file=paste("BP_down",getfiles[i],  sep=""))
  write.csv(ego6, file=paste("MF_down",getfiles[i],  sep=""))
  
  
}



############Venn diagrams############
# this analysis is based on differentially expressed genes between the reference (flat, unstimulated) and the other groups
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
library(readxl)
DE_US_porous_US_flat<-read.csv("DE_genes_porous_flat_US.csv")
DE_AB_flat_US_flat<-read.csv("DE_genes_AB_US_flat.csv")
DE_AB_porous_US_flat<-read.csv("DE_genes_AB_porous_US_flat.csv")


#upregulated genes
DE_US_porous_US_flat<-DE_US_porous_US_flat[DE_US_porous_US_flat$logFC>0, ]
DE_AB_flat_US_flat<-DE_AB_flat_US_flat[DE_AB_flat_US_flat$logFC>0, ]
DE_AB_porous_US_flat<-DE_AB_porous_US_flat[DE_AB_porous_US_flat$logFC>0, ]

set1<-DE_US_porous_US_flat$X
set2<-DE_AB_flat_US_flat$X
set3<-DE_AB_porous_US_flat$X

all_DE_genes<-Reduce(union, list(set1,set2,set3))
pheatmap(logCPM[all_DE_genes,], scale="row")
write.csv(all_DE_genes, file="all_DE_genes.csv", quote =FALSE)

get_genes_all3<-Reduce (intersect, list(set1,set2,set3))
get_genes_13<-setdiff(intersect(set1,set3),get_genes_all3)
get_genes_12<-setdiff(intersect(set1,set2),get_genes_all3)
get_genes_23<-setdiff(intersect(set2,set3),get_genes_all3)
get_genes_1<-setdiff(set1,union(set2,set3))
get_genes_2<-setdiff(set2,union(set1,set3))
get_genes_3<-setdiff(set3,union(set2,set1))

write.csv(get_genes_all3, file="get_genes_all3.csv", quote =FALSE)
write.csv(get_genes_13, file="get_genes_13.csv", quote =FALSE)
write.csv(get_genes_12, file="get_genes_12.csv", quote =FALSE)
write.csv(get_genes_23, file="get_genes_23.csv", quote =FALSE)
write.csv(get_genes_1, file="get_genes_1.csv", quote =FALSE)
write.csv(get_genes_2, file="get_genes_2.csv", quote =FALSE)
write.csv(get_genes_3, file="get_genes_3.csv", quote =FALSE)

x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_all3)
top10_all3_1<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_all3)
top10_all3_2<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_all3)
top10_all3_3<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_13)
top10_13_1<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_13)
top10_13_3<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_12)
top10_12_1<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_12)
top10_12_2<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_23)
top10_23_2<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_23)
top10_23_3<-x[order(x$PValue, decreasing = FALSE),]


x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_1)
top10_1<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_2)
top10_2<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_3)
top10_3<-x[order(x$PValue, decreasing = FALSE),]


get_top10genes<-Reduce(union, list(top10_1$X[1:10],top10_2$X[1:10],top10_3$X[1:10],
                                   top10_12_1$X[1:10],top10_12_2$X[1:10],
                                   top10_13_1$X[1:10],top10_13_3$X[1:10],
                                   top10_23_2$X[1:10],top10_23_3$X[1:10],
                                   top10_all3_1$X[1:10],top10_all3_2$X[1:10],top10_all3_3$X[1:10]))

write.csv(get_top10genes,file="get_top10genes.csv",quote = FALSE)
write.csv(top10_1,file="get_top10_1.csv",quote = FALSE)
write.csv(top10_2,file="get_top10_2.csv",quote = FALSE)
write.csv(top10_3,file="get_top10_3.csv",quote = FALSE)

write.csv(top10_12_1,file="get_top10_12_1.csv",quote = FALSE)
write.csv(top10_12_2,file="get_top10_12_2.csv",quote = FALSE)

write.csv(top10_23_2,file="get_top10_23_2.csv",quote = FALSE)
write.csv(top10_23_3,file="get_top10_23_3.csv",quote = FALSE)

write.csv(top10_13_1,file="get_top10_13_1.csv",quote = FALSE)
write.csv(top10_13_3,file="get_top10_13_3.csv",quote = FALSE)

write.csv(top10_all3_1,file="get_top10_all3_1.csv",quote = FALSE)
write.csv(top10_all3_2,file="get_top10_all3_2.csv",quote = FALSE)
write.csv(top10_all3_3,file="get_top10_all3_3.csv",quote = FALSE)

venn.diagram(
  x = list(
    DE_US_porous_US_flat$X , 
    DE_AB_flat_US_flat$X, 
    DE_AB_porous_US_flat$X
  ),
  category.names = c("US_porous VS US_flat" , "AB_flat VS US_flat" , "AB_porous VS US_flat"),
  #filename=VennDiagram,
  #filename= NULL,
  filename = 'venn_upregulated.png',
  output = TRUE,
  #imagetype="png" ,
  #height = 480 , 
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 1,
  fontfamily = "sans",
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)



#downregulated genes
DE_US_porous_US_flat<-read.csv("DE_genes_porous_flat_US.csv")
DE_AB_flat_US_flat<-read.csv("DE_genes_AB_US_flat.csv")
DE_AB_porous_US_flat<-read.csv("DE_genes_AB_porous_US_flat.csv")


DE_US_porous_US_flat<-DE_US_porous_US_flat[DE_US_porous_US_flat$logFC<0, ]
DE_AB_flat_US_flat<-DE_AB_flat_US_flat[DE_AB_flat_US_flat$logFC<0, ]
DE_AB_porous_US_flat<-DE_AB_porous_US_flat[DE_AB_porous_US_flat$logFC<0, ]


set1<-DE_US_porous_US_flat$X
set2<-DE_AB_flat_US_flat$X
set3<-DE_AB_porous_US_flat$X

all_DE_genes<-Reduce(union, list(set1,set2,set3))
pheatmap(logCPM[all_DE_genes,], scale="row")
write.csv(all_DE_genes, file="all_DE_genes.csv", quote =FALSE)

get_genes_all3<-Reduce (intersect, list(set1,set2,set3))
get_genes_13<-setdiff(intersect(set1,set3),get_genes_all3)
get_genes_12<-setdiff(intersect(set1,set2),get_genes_all3)
get_genes_23<-setdiff(intersect(set2,set3),get_genes_all3)
get_genes_1<-setdiff(set1,union(set2,set3))
get_genes_2<-setdiff(set2,union(set1,set3))
get_genes_3<-setdiff(set3,union(set2,set1))


write.csv(get_genes_all3, file="get_genes_all3_down.csv", quote =FALSE)
write.csv(get_genes_13, file="get_genes_13_down.csv", quote =FALSE)
write.csv(get_genes_12, file="get_genes_12_down.csv", quote =FALSE)
write.csv(get_genes_23, file="get_genes_23_down.csv", quote =FALSE)
write.csv(get_genes_1, file="get_genes_1_down.csv", quote =FALSE)
write.csv(get_genes_2, file="get_genes_2_down.csv", quote =FALSE)
write.csv(get_genes_3, file="get_genes_3_down.csv", quote =FALSE)


x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_all3)
top10_all3_1<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_all3)
top10_all3_2<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_all3)
top10_all3_3<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_13)
top10_13_1<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_13)
top10_13_3<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_12)
top10_12_1<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_12)
top10_12_2<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_23)
top10_23_2<-x[order(x$PValue, decreasing = FALSE),]
x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_23)
top10_23_3<-x[order(x$PValue, decreasing = FALSE),]


x<-subset(DE_US_porous_US_flat, DE_US_porous_US_flat$X %in% get_genes_1)
top10_1<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_AB_flat_US_flat, DE_AB_flat_US_flat$X %in% get_genes_2)
top10_2<-x[order(x$PValue, decreasing = FALSE),]

x<-subset(DE_AB_porous_US_flat, DE_AB_porous_US_flat$X %in% get_genes_3)
top10_3<-x[order(x$PValue, decreasing = FALSE),]


get_top10genes<-Reduce(union, list(top10_1$X[1:10],top10_2$X[1:10],top10_3$X[1:10],
                                   top10_12_1$X[1:10],top10_12_2$X[1:10],
                                   top10_13_1$X[1:10],top10_13_3$X[1:10],
                                   top10_23_2$X[1:10],top10_23_3$X[1:10],
                                   top10_all3_1$X[1:10],top10_all3_2$X[1:10],top10_all3_3$X[1:10]))

write.csv(get_top10genes,file="get_top10genes_down.csv",quote = FALSE)
write.csv(top10_1,file="get_top10_1_down.csv",quote = FALSE)
write.csv(top10_2,file="get_top10_2_down.csv",quote = FALSE)
write.csv(top10_3,file="get_top10_3_down.csv",quote = FALSE)

write.csv(top10_12_1,file="get_top10_12_1_down.csv",quote = FALSE)
write.csv(top10_12_2,file="get_top10_12_2_down.csv",quote = FALSE)

write.csv(top10_23_2,file="get_top10_23_2_down.csv",quote = FALSE)
write.csv(top10_23_3,file="get_top10_23_3_down.csv",quote = FALSE)

write.csv(top10_13_1,file="get_top10_13_1_down.csv",quote = FALSE)
write.csv(top10_13_3,file="get_top10_13_3_down.csv",quote = FALSE)

write.csv(top10_all3_1,file="get_top10_all3_1_down.csv",quote = FALSE)
write.csv(top10_all3_2,file="get_top10_all3_2_down.csv",quote = FALSE)
write.csv(top10_all3_3,file="get_top10_all3_3_down.csv",quote = FALSE)



venn.diagram(
  x = list(
    DE_US_porous_US_flat$X , 
    DE_AB_flat_US_flat$X, 
    DE_AB_porous_US_flat$X
  ),
  category.names = c("US_porous VS US_flat" , "AB_flat VS US_flat" , "AB_porous VS US_flat"),
  #filename=VennDiagram,
  #filename= NULL,
  filename = 'venn_downregulated.png',
  output = TRUE,
  #imagetype="png" ,
  #height = 480 , 
  #width = 480 , 
  #resolution = 300,
  #compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 1,
  fontfamily = "sans",
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)





############Perform GO enrichment on the different gene subsets from Venn diagrams (figures 1 and S1).############ 
#The gene lists are in the folder "Venn_genes".

getfiles<-list.files("Venn_genes")
for (i in 1:length(getfiles)){
  DE_genes<-read.csv(paste("Venn_genes/",getfiles[i],sep=""))  
  
  
  
  gene.df <- bitr(DE_genes$x, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  
  ego <- enrichGO(gene         = gene.df$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  # keytype       = 'SYMBOL',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable =TRUE)

  ego2 <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   # keytype       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable = TRUE)
  
  
  ego3 <- enrichGO(gene         = gene.df$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   # keytype       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable=TRUE)
  
  
  
  
  pdf(paste("Path_en_",getfiles[i],"_Venn.pdf",sep=""))
  print(dotplot(ego,showCategory=20))
  print(dotplot(ego2,showCategory=20))
  print(dotplot(ego3,showCategory=20))
  
  
  
  
  write.csv(ego, file=paste("CC_",getfiles[i], sep=""))
  write.csv(ego2, file=paste("BP_",getfiles[i], sep=""))
  write.csv(ego3, file=paste("MF_",getfiles[i], sep=""))
  
  
  dev.off()
  
}



############Compare stimulated cells on porous VS flat surface############
get_genes_Ab_surf<-read.csv("DE_genes_porous_flat_AB.csv")
ordered_genes<-get_genes_Ab_surf$logFC
names(ordered_genes)<-get_genes_Ab_surf$X


#Order genes based on logFC between porous and flat stimulated
ordered_genes<-topTags(res_porous_flat_AB,n=Inf)[,1]$table
test<-ordered_genes$logFC
test<-2^test
names(test)<-rownames(ordered_genes)
test<-test[order(test,decreasing = TRUE)]


library(fgsea)
library("GSEABase")
library(qusage)
#read gmt file(s) downloaded from MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/index.jsp

gmtfile <- "c2.all.v7.0.symbols.gmt"

# other tested gmt files:  h.all.v6.1.symbols.gmt, c5.all.v7.0.symbols.gmt, c7.all.v7.0.symbols.gmt


react_pathways <- read.gmt(gmtfile)

fgseaRes <- fgsea(pathways = react_pathways, 
                  stats = test,
                  minSize=5,
                  maxSize=2000,
                  nperm=10000)
head(fgseaRes[order(pval), ])
test2<-fgseaRes
test2$leadingEdge<-as.character(test2$leadingEdge)


write.csv(test2, file="fgsea_c2.csv", quote=FALSE)


#plot MAPK pathway

pathway<-"REACTOME_MAPK_FAMILY_SIGNALING_CASCADES"

pdf(file=paste("MAPK_porous_pos_VS_flat_pos.pdf",sep=""),height = 30, width=56)
plotEnrichment(react_pathways[[pathway]],test, ticksSize = 1) +
  labs(title=pathway)
dev.off()




############Compare unstimulated cells on porous VS flat surface############
get_genes_US_surf<-read.csv("DE_genes_porous_flat_US.csv")
ordered_genes<-get_genes_US_surf$logFC
names(ordered_genes)<-get_genes_US_surf$X


ordered_genes<-topTags(res_porous_flat_US,n=Inf)[,1]$table
test<-ordered_genes$logFC
test<-2^test
names(test)<-rownames(ordered_genes)
test<-test[order(test,decreasing = TRUE)]


library(fgsea)
library("GSEABase")
library(qusage)

gmtfile <- "c2.all.v7.0.symbols.gmt"
react_pathways <- read.gmt(gmtfile)

fgseaRes <- fgsea(pathways = react_pathways, 
                  stats = test,
                  minSize=5,
                  maxSize=2000,
                  nperm=10000)
head(fgseaRes[order(pval), ])
test2<-fgseaRes
test2$leadingEdge<-as.character(test2$leadingEdge)


write.csv(test2, file="fgsea_c2.csv", quote=FALSE)



pathway<-"REACTOME_MAPK_FAMILY_SIGNALING_CASCADES"
pdf(file=paste("MAPK_porous_neg_VS_flat_neg.pdf",sep=""),height = 30, width=56)
plotEnrichment(react_pathways[[pathway]],test, ticksSize = 1) +
  labs(title=pathway)
dev.off()


#### Hallmark genes for MAPK pathway 
mapk_genes<-unique(c("NRG4", "CSF2", "IL2", "DUSP6", "IL3", "DUSP4", "ACTN2", "AREG", "DUSP5", "PDGFA", "IL2RA", 
                     "MYC", "DUSP8", "DUSP2", "JAK2", "SPTBN4", "PEA15", "KSR1", "CDC42EP3", "AGO2", "PSME3", 
                     "PSMD14", "SPRED1", "PHB", "MAPK6", "PSMA6", "PSMB2", "PSMB5", "CAMK2D", "PSMC4", "PSMA3", 
                     "PSMD12", "APBB1IP", "PSMA1", "CNKSR1", "FYN", "PSMA2", "PSMD1", "PSMB7", "PPP5C", "PSMD11", 
                     "HBEGF", "PSMB3", "DUSP10", "PSMD8", "MAPKAPK5", "PTPN11", "FGFR4", "PSMD9", "PSMC1", "DUSP1", 
                     "PSMC3", "PSMA7", "PSMB6", "JUN",
                     "MAP2K4","MAP3K14","MAP3K8","MAPKAPK3", "MAP2K1","MAP3K14", "MAP3K8","MAP2K3"
                    
))


mapk_genes<-intersect(mapk_genes,rownames(logCPM))


Tcell_activ<-unique(c("CD45RO","PTPRC","PTPRC","CD44","CD69","IL2RA","NR4A1","CD8A","IRF4","NFKBIB","NFKBIA",
                      "JUN","NFATC1","CD7","LILR4B","CD200R1","PRDM1","CD164","CD160","CD244",
                      "LAT","CD27","FOSL2","JAK1","TBX21","EOMES", "TCF7","ZEB2","EZH2",
                      "FOS","CD3","PDCD1","CTLA4","HAVCR2","IFNG","TNFA","IL2","KLRG1","SELL","CX3CR1","GPR13",
                      "TIGIT","ENTPD1","CCR7","TOX","CD45","SELL","IL7R","HLA-DRA","HLA-A","HLA-E", "NFATC1", "LAT",
                      "IL2","CD69", "PTPN22","RELA","MAP3K8", "RELA","TRAV8-4","TRAV29DV5","TRBV7-9","CD3D","CD3E",
                      "CALM2","BCL2"))

Tcell_activ<-intersect(Tcell_activ,rownames(logCPM))

col_tags<-as.data.frame(colnames(logCPM))
col_tags$stimulation<-c(rep("-",6),rep("+",6))
col_tags$surface<-c(rep("flat",3),rep("porous",6),rep("flat",3))
rownames(col_tags)<-col_tags$`colnames(logCPM)`
col_tags<-col_tags[,2:3]


my_colour = list(
  surface = c(flat = "darkgrey", porous = "black"),
  stimulation = c('-' = "darkcyan", '+' = "darkorchid4")
)



pdf(file=paste("MAPK_heatmap_all.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[mapk_genes,],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()



pdf(file=paste("Tcellactiv_selected.pdf",sep=""),height = 15, width=8)
pheatmap(logCPM[Tcell_activ,],scale = "row", fontsize_row = 10, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()


############Hypergeometric plots (figure 1)############
#read list of enriched pathways from file
stim<-read_xlsx ("enriched_pathways.xlsx")
stim<-as.data.frame(stim)
stim<-subset(stim, stim$padj<0.05)
stim<-stim[order(stim$size),]

stim$pathway<-factor(stim$pathway, levels=stim$pathway)

library(extrafont)
loadfonts(device = "win")
fonts$register_font('Arial')


p<-ggplot()+
  geom_point(data=stim, mapping=aes(x=NES, y=pathway, size=size, colour = -log10(padj))) +
  scale_color_gradient2(midpoint = mean(-log10(stim$padj)),low="blue", mid="orange",high="red" )+

  theme_bw()+

  scale_size(range = c(1, 3))+
  theme(text=element_text(size=4,  family="Arial"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggsave(p, filename = "gsea_por_pos_VS_flat_pos_padj05_all_new.pdf",height=4.8, width=7.8,
       units="cm", device = cairo_pdf)




#read list of enriched goterms from file

library(readxl)
library(ggplot2)
library(ggrepel)
stim<-read_xlsx ("GO_term_porous_flat_pos_sel.xlsx")
stim<-as.data.frame(stim)
stim$`up/down`<-as.factor(stim$`up/down`)
stim<-stim[order(stim$`ratio pathway`),]
stim$Description<-factor(stim$Description, levels=stim$Description)

pdf(file="go_term_por_pos_VS_flat_pos_sel.pdf",height = 8, width=10)
ggplot()+
  geom_point(data=stim, mapping=aes(x=-log10(qvalue), y=Description, size=`ratio pathway`, colour= `up/down`)) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#general heatmap - read differentially expressed genes
get_genes<-read.csv("DE_genes_porous_flat_AB.csv")
get_genes<-get_genes[order(get_genes$PValue),]

#load object y (edgeR object containing normalized and fitted data)
load("y.Robj")
logCPM <- cpm(y, prior.count=10, log=TRUE)
col_tags<-as.data.frame(colnames(logCPM))
col_tags$stimulation<-c(rep("-",6),rep("+",6))
col_tags$surface<-c(rep("flat",3),rep("porous",6),rep("flat",3))
rownames(col_tags)<-col_tags$`colnames(logCPM)`
col_tags<-col_tags[,2:3]


my_colour = list(
  surface = c(flat = "darkgrey", porous = "black"),
  stimulation = c('-' = "darkcyan", '+' = "darkorchid4")
)
library(pheatmap)

pdf(file=paste("general_heatmap_top200.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[as.character(get_genes$X)[1:200],],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

write.csv(get_genes[1:200,], file="top_200_genes.csv", quote=FALSE)


############MAPK pathways heatmaps (figure 3)############

library(edgeR)

#load edgeR object (y) containing preprocessed and fitted data, the fit, and the design matrix
load("y.Robj")

load("design.Robj")
load("fit.Robj")


library(fgsea)
library("GSEABase")
library(qusage)

gmtfile <- "c2.all.v7.0.symbols.gmt"

react_pathways <- read.gmt(gmtfile)

get_genes_MAPK<-react_pathways$REACTOME_MAPK_FAMILY_SIGNALING_CASCADES
get_genes_MAPK_targets<-react_pathways$YOSHIMURA_MAPK8_TARGETS_UP

logCPM <- cpm(y, prior.count=10, log=TRUE)
col_tags<-as.data.frame(colnames(logCPM))
col_tags$stimulation<-c(rep("-",6),rep("+",6))
col_tags$surface<-c(rep("flat",3),rep("porous",6),rep("flat",3))
rownames(col_tags)<-col_tags$`colnames(logCPM)`
col_tags<-col_tags[,2:3]


my_colour = list(
  surface = c(flat = "darkgrey", porous = "black"),
  stimulation = c('-' = "darkcyan", '+' = "darkorchid4")
)
library(pheatmap)

logCPM_noclust<-logCPM[,c(1:6,10:12,7:9)]

get_genes_MAPK<-intersect(get_genes_MAPK,rownames(logCPM))
get_genes_MAPK_targets<-intersect(get_genes_MAPK_targets,rownames(logCPM))


pdf(file=paste("REACTOME_MAPK_FAMILY_SIGNALING_CASCADES_all_clust.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[get_genes_MAPK,],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

pdf(file=paste("REACTOME_MAPK_FAMILY_SIGNALING_CASCADES_all_NOclust.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM_noclust[get_genes_MAPK,],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         cluster_cols = FALSE,
         annotation_colors = my_colour)
dev.off()




pdf(file=paste("YOSHIMURA_MAPK8_TARGETS_UP_all_clust.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[get_genes_MAPK_targets,],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

pdf(file=paste("YOSHIMURA_MAPK8_TARGETS_UP_all_NOclust.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM_noclust[get_genes_MAPK_targets,],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         cluster_cols = FALSE,
         annotation_colors = my_colour)
dev.off()

#Show MAPK genes which are differentially expressed between porous (stimulated) and flat (stimulated) conditions
get_genes_Ab_surf<-read.csv("DE_genes_porous_flat_AB.csv")


get_top_50<-subset(get_genes_Ab_surf,get_genes_Ab_surf$X %in% get_genes_MAPK)
get_top_50<-subset(get_top_50, get_top_50$logFC>0.5)
get_top_50<-get_top_50[order(get_top_50$PValue),]
pdf(file=paste("REACTOME_MAPK_FAMILY_SIGNALING_CASCADES_top50_UP_clust.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[as.character(get_top_50$X[1:50]),],scale = "row", fontsize_row = 5, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

write.csv(get_top_50, file="REACTOME_MAPK_FAMILY_SIGNALING_CASCADES.csv", quote=FALSE)



get_top_50<-subset(get_genes_Ab_surf,get_genes_Ab_surf$X %in% get_genes_MAPK_targets)
get_top_50<-subset(get_top_50, get_top_50$logFC>0.5)
get_top_50<-get_top_50[order(get_top_50$PValue),]
pdf(file=paste("YOSHIMURA_MAPK8_TARGETS_top50_UP_clust.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[as.character(get_top_50$X[1:50]),],scale = "row", fontsize_row = 5, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

write.csv(get_top_50, file="YOSHIMURA_MAPK8_TARGETS.csv", quote=FALSE)







############hypergeometric plot (figure 4)############

library(readxl)
library(ggplot2)
#read enriched pathways (from GSEA analysis) in porous positive versus flat positive 
stim<-read_xlsx ("enriched_pathways.xlsx")
stim<-as.data.frame(stim)
stim<-subset(stim, stim$padj<0.05)
stim<-stim[order(stim$size),]

stim$pathway<-factor(stim$pathway, levels=stim$pathway)


stim$pathway<-factor(stim$pathway, levels=stim$pathway)

pdf(file="gsea_por_pos_VS_flat_pos_padj.pdf",height = 1.88976, width=2.95276) 
ggplot()+
  scale_size(range = c(1, 3))+
  geom_point(data=stim, mapping=aes(x=NES, y=pathway, size=size, colour = -log10(padj))) +
  scale_color_gradient2(midpoint = mean(-log10(stim$padj)),low="blue", mid="orange",high="red")+
  theme_bw()+
  xlim(1.2,1.8)+
  theme(text = element_text(size=5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.1),
        legend.key.size = unit(0.3, "cm")
        
  )

dev.off()

#read file containing enriched GO terms in porous negative VS flat negative
stim<-read_xlsx ("GO_term_porous_flat_US.xlsx")
stim<-as.data.frame(stim)
stim$`up/down`<-as.factor(stim$`up/down`)
stim<-stim[order(stim$`ratio pathway`),]
stim$Description<-factor(stim$Description, levels=stim$Description)

pdf(file="go_term_por_neg_VS_flat_neg_sel2.pdf",height = 8, width=10)
ggplot()+
  geom_point(data=stim, mapping=aes(x=-log10(qvalue), y=Description, size=`ratio pathway`, colour= `up/down`)) +
  scale_colour_brewer(palette = "Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()





############TCR signaling heatmap (figure 4)############

library(fgsea)
library("GSEABase")
library(qusage)
library(pheatmap)

gmtfile <- "c5.all.v7.0.symbols.gmt"

react_pathways <- read.gmt(gmtfile)


get_genes_TCR<-react_pathways$GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY

#read differentially expressed genes between porous negative and flat negative
get_genes_US_surf<-read.csv("DE_genes_porous_flat_US.csv")

#load edgeR object (y) containing normalized and fitted data
load("y.Robj")
logCPM <- cpm(y, prior.count=10, log=TRUE)
col_tags<-as.data.frame(colnames(logCPM))
col_tags$stimulation<-c(rep("-",6),rep("+",6))
col_tags$surface<-c(rep("flat",3),rep("porous",6),rep("flat",3))
rownames(col_tags)<-col_tags$`colnames(logCPM)`
col_tags<-col_tags[,2:3]


my_colour = list(
  surface = c(flat = "darkgrey", porous = "black"),
  stimulation = c('-' = "darkcyan", '+' = "darkorchid4")
)

#select upregulated genes from the signaling pathway
get_genes_TCR<-intersect(get_genes_TCR,rownames(logCPM))


pdf(file=paste("GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY_all.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[get_genes_TCR,1:6],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

#plot top 50 genes (ordered by pvalue)
get_top_50<-subset(get_genes_US_surf,get_genes_US_surf$X %in% get_genes_TCR)
get_top_50<-subset(get_top_50, get_top_50$logFC>0.5)
get_top_50<-get_top_50[order(get_top_50$PValue),]
pdf(file=paste("GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY_top50_UP_clust.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[as.character(get_top_50$X[1:50]),1:6],scale = "row", fontsize_row = 7, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

write.csv(get_top_50, file="T_CELL_RECEPTOR_SIGNALING_PATHWAY.csv", quote=FALSE)


############Heatmap of top 200 differentially expressed genes between porous negative and flat negative############
library(pheatmap)
load("edgeR")
load("ggplot2")
#load edgeR object (y) containing normalized and fitted data
load("y.Robj")

get_genes_US_surf<-read.csv("DE_genes_porous_flat_US.csv")
logCPM <- cpm(y, prior.count=10, log=TRUE)
col_tags<-as.data.frame(colnames(logCPM))
col_tags$stimulation<-c(rep("-",6),rep("+",6))
col_tags$surface<-c(rep("flat",3),rep("porous",6),rep("flat",3))
rownames(col_tags)<-col_tags$`colnames(logCPM)`
col_tags<-col_tags[,2:3]


my_colour = list(
  surface = c(flat = "darkgrey", porous = "black"),
  stimulation = c('-' = "darkcyan", '+' = "darkorchid4")
)


pdf(file=paste("general_heatmap_top200.pdf",sep=""),height = 21, width=8)
pheatmap(logCPM[as.character(get_genes$X)[1:200],],scale = "row", fontsize_row = 3, annotation_col = col_tags,
         annotation_colors = my_colour)
dev.off()

write.csv(get_genes[1:200,], file="top_200_genes.csv", quote=FALSE)



############3D plots############
library(edgeR)
library(rgl)
library("plot3D")
#load edgeR object (y) with normalized and fitted data, design matrix, and fitting results 
load("y.Robj")

load("design.Robj")
load("fit.Robj")

get_dim<-p@.Data[[3]]

mds.1<-get_dim[,1]
mds.2<-get_dim[,2]
mds.3<-get_dim[,3]

pdf(file="MDS_3D.pdf",height = 7, width=10)
scatter3D(mds.1, mds.2, mds.3,
          colvar = NULL, bty = "g", pch=20, cex=5,
          xlab = c("MDS 1"), ylab="MDS 2",zlab ="MDS 3",
          #phi=30, theta =30,
          #as.factor(c(rep("flat negative",3), rep("porous negative",3),rep("porous positive",3), rep("flat positive",3))),
          col = c(rep("#1B9E77",3), rep("#D95F02",3), rep("#7570B3",3),rep("black",3)))


legend("topright", c("flat negative","porous negative", "porous positive", "flat positive"), pch = 16, 
       col = c("#1B9E77","#D95F02","#7570B3","black"), cex=1, inset=c(0.00000005,0.7))

dev.off()
############Volcano plots############
library(EnhancedVolcano)
#load R object containing differentially expressed genes between porous positive and flat positive (columns contain logFC and adjusted pvalue (FDR))
load("volcano_porous_flat_AB.Robj")
pdf(file="volcano_porous_pos_VS_flat_pos.pdf",height = 8, width=10)
EnhancedVolcano(volcano_data,
                lab = rownames(volcano_data),
                ylim = c(0,10),
                selectLab = c("IFNG", "CD69","PTPRC","CD25","TUBA","IL2RA","TUBA1A","TUBA1B","IL2",
                              "NR4A1","ITGAL","NFKBIZ","NFKBIA","TUBB","PTPRCAP"),
                subtitle = 'Porous + VS flat +',
                FCcutoff = 0.5,
                x = 'logFC',
                y = 'FDR',
                transcriptLabSize = 4.0,
                pCutoff =5*10e-2,
                # xlim = c(-5, 8)
)
dev.off()

