#聚类
data<-read.table("uniSigExp.txt",head=T,sep='\t',check.names=F,row.names=1)
ata<-as.matrix(data)
maxK=9
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")

Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)}#end for i# The optimal K
optK = Kvec[which.min(PAC)]#输出最佳K值
#输出结果
clusterNum=optK                  #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="cluster.txt",sep="\t",quote=F,col.names=F)

#Survival
library(survival)
library("survminer")
rt=read.table("cluster.txt",header=T,sep="\t");rt$futime=rt$futime/12
diff=survdiff(Surv(futime,fustat) ~Gene_cluster,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime,fustat) ~ Gene_cluster, data = rt)
ggsurvplot(fit, 
           data=rt,
           conf.int=T,
           pval=paste0 ("P = ",pValue),
           pval.size=6,
           legend.title="Gene_cluster",
           legend.labs=levels(factor(rt[,"Gene_cluster"])),
           legend = c(0.9, 0.9),
           font.legend=10,
           xlab="Time(years)",
           break.time.by = 1,
           palette = "nejm",
           surv.median.line = "hv",
           risk.table=T,conf.int.style = "step",
           cumevents=F,ncensor.plot = TRUE,ncensor.plot.height = 0.25,
           risk.table.height=.30)
		   
		   
#Immune
library(utils)
library(GSVA)
library(ComplexHeatmap) 
source("pheatmap_translate.R")
library(circlize)
library(viridis)
library(gplots)
library(data.table)
library(estimate)
source("annTrackScale.R") 
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 


standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
    outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
    if (!is.null(halfwidth)) {
        outdata[outdata>halfwidth]=halfwidth
        outdata[outdata<(-halfwidth)]= -halfwidth
    }
    return(outdata)
}

annCol.tcga<-read.table("Cluster.txt",head=T,sep='\t',check.names = F,row.names = 1)
tpm <-read.table("E-MATB-1980.txt",head=T,sep='\t',check.names = F,row.names = 1)
tpm<-tpm[,intersect(row.names(annCol.tcga),colnames(tpm))]
immune.signature <- read.table("Curated_Immune_Cell_Signature.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)


cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
    immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}


imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4") 


immune.sig.ccr.order <- c("T.cells.CD8", 
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")




indata<-tpm
write.table(indata,file = "TCGA_log2TPM_hugo.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

filterCommonGenes(input.f = "TCGA_log2TPM_hugo.txt" , output.f="TCGA_log2TPM_hugo_ESTIMATE.txt", id="GeneSymbol")

estimateScore("TCGA_log2TPM_hugo_ESTIMATE.txt","TCGA_log2TPM_hugo_estimate_score.txt", platform="affymetrix")

est.tcga <- read.table(file = "TCGA_log2TPM_hugo_estimate_score.txt",header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.tcga) <- est.tcga[,2]; colnames(est.tcga) <- est.tcga[1,]; est.tcga <- est.tcga[-1,c(-1,-2)];
est.tcga <- sapply(est.tcga, as.numeric); rownames(est.tcga) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.tcga.backup = as.data.frame(est.tcga); colnames(est.tcga.backup) <- colnames(indata)


est.tcga <- annTrackScale(indata = est.tcga, halfwidth = 2, poolsd = F); est.tcga <- as.data.frame(t(est.tcga)) 
rownames(est.tcga) <- colnames(tpm)

tcga.immune.gsva <- gsva(as.matrix(tpm),
                         immune.sig.ccr,
                         method = "gsva")





clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
blue <- "#5bc0eb"
gold <- "#ECE700"
cyan <- "#00B3D0"

annCol.tcga$ImmuneScore <- as.numeric(est.tcga[rownames(annCol.tcga),"ImmuneScore"])
annCol.tcga$StromalScore <- as.numeric(est.tcga[rownames(annCol.tcga),"StromalScore"])
annCol.tcga <- annCol.tcga[order(annCol.tcga$Cluster),]
annColors.tcga <- list() 
annColors.tcga[["Cluster"]] <- c("IR1" = clust.col[1],
                               "IR2" = clust.col[2],
                               "IR3" = clust.col[3]
                               )
annColors.tcga[["ImmuneScore"]] <- inferno(64)
annColors.tcga[["StromalScore"]] <- viridis(64)


indata <- tpm[intersect(rownames(tpm),imm.targets),]
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.tcga)],halfwidth = 2), 
                border_color = NA, 
                annotation_col = annCol.tcga[,c("Cluster","StromalScore","ImmuneScore")],
                annotation_colors = annColors.tcga[c("Cluster","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T, 
                show_colnames = F, 
                cellheight = 12, 
                cellwidth = 0.6, 
                name = "ICI", 
                cluster_rows = F, 
                cluster_cols = F)


hm2<-pheatmap(standarize.fun(tcga.immune.gsva[immune.sig.ccr.order,rownames(annCol.tcga)],halfwidth = 1), 
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22), 
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "TIME",
                cluster_rows = F,
                cluster_cols = F)
draw(hm1 %v% hm2, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")