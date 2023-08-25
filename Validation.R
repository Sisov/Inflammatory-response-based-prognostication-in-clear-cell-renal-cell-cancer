library("MOVICS")

dat<-read.table("uniq.symbol.txt",head=T,sep='\t',check.names = F,row.names = 1)
survinfo<-read.table("cluster.txt",head=T,sep='\t',check.names = F,row.names = 1)

identical(colnames(dat),row.names(survinfo))

pseudo.moic.res<- list("clust.res" = survinfo,
                                        "mo.method" = "Consensus") 
fpkm<-dat
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)

runDEA(dea.method = "limma",
       expr       = fpkm, 
       moic.res   = pseudo.moic.res,
       prefix     = "TCGA-ccRCC")
	   
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-ccRCC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = fpkm, # use normalized expression as heatmap input
#                       annCol        = annCol, 
#                       annColors     = annColors, 
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")


marker.down <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-ccRCC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "down", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = fpkm, # use normalized expression as heatmap input
#                       annCol        = annCol, # 
#                       annColors     = annColors, # 
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")


MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
#up regulation biomarker
gsea.up <- runGSEA(moic.res     = pseudo.moic.res,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "TCGA-ccRCC", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")

#down regulation biomarker
gsea.dn <- runGSEA(moic.res     = pseudo.moic.res,
                   dea.method   = "limma",
                   prefix       = "TCGA-ccRCC",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP") 

mRNA<-read.table("E-MATB-1980.txt",head=T,sep='\t',check.names = F,row.names = 1)
surv.info<-read.table("E-MTAB.txt",head=T,sep='\t',check.names = F,row.names = 1)
identical(row.names(surv.info),colnames(mRNA))
validation<-list("validation" = mRNA,"clin.info"=surv.info)
#NTP Evaluation
yau.ntp.pred <- runNTP(expr       = validation$validation,
                       templates  = marker.up$templates, 
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU") 

table(pseudo.moic.res$clust.res$clust)
#PAM Evaluation
yau.pam.pred <- runPAM(train.expr  = fpkm,
                       moic.res    = pseudo.moic.res,
                       test.expr   = validation$validation)

print(yau.pam.pred$IGP)

tcga.ntp.pred <- runNTP(expr      = fpkm,
                        templates = marker.up$templates,
                        doPlot    = FALSE) 
						
tcga.pam.pred <- runPAM(train.expr  = fpkm,
                        moic.res    = pseudo.moic.res,
                        test.expr   = fpkm)


runKappa(subt1     = pseudo.moic.res$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")

runKappa(subt1     = pseudo.moic.res$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")

runKappa(subt1     = yau.pam.pred$clust.res$clust
         subt2     = yau.ntp.pred$clust.res$clust
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR YAU")