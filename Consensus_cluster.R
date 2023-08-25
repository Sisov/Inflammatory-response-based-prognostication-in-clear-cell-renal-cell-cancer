#Clustering
gene<-read.table("gene.txt",header=F,sep="\t",check.names=F)
dat<-read.table("uniq.symbol.txt",head=T,sep='\t',check.names=F,row.names=1)
expr<-dat[intersect(gene$VA,row.names(dat),]
time<-read.table("time.txt",head=T,sep='\t',check.names=F,row.names=1)
expr<-expr[,intersect(row.names(time),colnames(expr)]

library(survival)
pFilter=0.05                                                    
rt=cbind(time,t(expr))    
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 if(coxP<pFilter){
     sigGenes=c(sigGenes,i)
		 outTab=rbind(outTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
  }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

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
clusterNum=optK                 
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

annTrackScale <- function(indata=NULL, halfwidth=NULL, poolsd=FALSE) {
  if ( !(is.vector(indata) | is.matrix(indata) | is.data.frame(indata)) ) {
    stop("indata must be vector, matrix or data frame!")
  }
  
  if (is.vector(indata)) {
    if (sd(indata)!=0) {
      outdata = indata-median(indata, na.rm=T)  
      
      if (poolsd==FALSE) {      
        pos <- outdata[outdata>=0]
        neg <- outdata[outdata<0]
        sdpos <- sd(pos, na.rm = T)
        sdneg <- sd(neg, na.rm = T)
        pos <- pos/sdpos
        neg <- neg/sdneg
        outdata <- c(pos, neg)[names(outdata)]
      }else{    
        std <- sd(outdata, na.rm = T)
        outdata <- outdata/std
      }
      
      if (!is.null(halfwidth)) {
        outdata[outdata>halfwidth]=halfwidth
        outdata[outdata<(-halfwidth)]= -halfwidth
        outdata[outdata<0] <- halfwidth*outdata[outdata<0]/abs(min(outdata))
        outdata[outdata>0] <- halfwidth*outdata[outdata>0]/max(outdata)
      }
    }else{
      outdata <- rep(0, times=length(indata))
    }
    
  }
  
  if (is.matrix(indata) | is.data.frame(indata)) {    
    outdata = sweep(indata,1, apply(indata,1,median,na.rm=T))
    for (m in 1:nrow(outdata)) {
      tmp <- as.numeric(outdata[m, ]); names(tmp) <- colnames(outdata)
      
      if(sd(tmp)!=0) {
        if (poolsd==FALSE) {
          pos <- tmp[tmp>=0]
          neg <- tmp[tmp<0]
          sdpos <- sd(pos, na.rm = T)
          sdneg <- sd(neg, na.rm = T)
          pos <- pos/sdpos
          neg <- neg/sdneg
          outdata[m,] <- c(pos, neg)[names(tmp)]       
        }else{
          std <- sd(tmp, na.rm = T)
          outdata[m,] <- tmp/std
        }      
      }else{
        outdata[m,] <- rep(0, times=length(tmp))
      }
    }    
    
    if (!is.null(halfwidth)) {
      outdata[outdata>halfwidth]=halfwidth
      outdata[outdata<(-halfwidth)]= -halfwidth
      for (k in 1:nrow(outdata)) {
        rowk <- outdata[k, ]
        rowk[rowk<0] <- halfwidth*rowk[rowk<0]/abs(min(rowk))
        rowk[rowk>0] <- halfwidth*rowk[rowk>0]/max(rowk)
        outdata[k,] <- rowk
      }
    }
  }
  
  return(outdata)
}
pheatmap = function(mat, 
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
    kmeans_k = NA, 
    breaks = NA, 
    border_color = ifelse(nrow(mat) < 100 & ncol(mat) < 100, "grey60", NA),
    cellwidth = NA, 
    cellheight = NA, 
    scale = "none", 
    cluster_rows = TRUE,
    cluster_cols = TRUE, 
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean", 
    clustering_method = "complete",
    clustering_callback = NA, 
    cutree_rows = NA, 
    cutree_cols = NA,
    treeheight_row = ifelse(class(cluster_rows) == "hclust" || cluster_rows, 50, 0), 
    treeheight_col = ifelse(class(cluster_cols) == "hclust" || cluster_cols, 50, 0), 
    legend = TRUE, 
    legend_breaks = NA,
    legend_labels = NA, 
    annotation_row = NA, 
    annotation_col = NA,
    annotation = NA, 
    annotation_colors = NA, 
    annotation_legend = TRUE,
    annotation_names_row = TRUE, 
    annotation_names_col = TRUE,
    drop_levels = TRUE, 
    show_rownames = TRUE, 
    show_colnames = TRUE, 
    main = NA,
    fontsize = 10, 
    fontsize_row = fontsize, 
    fontsize_col = fontsize,
    angle_col = c("270", "0", "45", "90", "315"), 
    display_numbers = FALSE,
    number_format = "%.2f", 
    number_color = "grey30", 
    fontsize_number = 0.8 * fontsize, 
    gaps_row = NULL, 
    gaps_col = NULL, 
    labels_row = NULL,
    labels_col = NULL, 
    filename = NA, 
    width = NA, 
    height = NA,
    silent = FALSE, 
    na_col = "#DDDDDD", 
    name = NULL,

    # argument specific for Heatmap()
    heatmap_legend_param = list(),
    ...,
    run_draw = FALSE
) {
  
    if(is.data.frame(mat)) {
        warning_wrap("The input is a data frame, convert it to the matrix.")
        mat = as.matrix(mat)
    }

    if(!identical(kmeans_k, NA)) {
        warning_wrap("argument `kmeans_k` is not suggested to use in pheatmap -> Heatmap translation because it changes the input matrix. You might check `row_km` and `column_km` arguments in Heatmap().")
        km = kmeans(mat, centers = kmeans_k)
        mat = km$centers
        rownames(mat) = paste0("Cluster: ", seq_along(km$size), ", Size: ", km$size)
    }

    if("row" %in% scale) {

        if(any(is.na(mat))) {
            mat = (mat - rowMeans(mat, na.rm = TRUE))/rowSds(mat, na.rm = TRUE)
        } else {
            mat = t(scale(t(mat)))
        }
    } else if("column" %in% scale) {
        if(any(is.na(mat))) {
            mat = t((t(mat) - colMeans(mat, na.rm = TRUE))/colSds(mat, na.rm = TRUE))
        } else {
            mat = scale(mat)
        }
    }

    ht_param = list(matrix = mat)

    if(!identical(scale, "none") && !identical(breaks, NA)) {
        warning_wrap("It not suggested to both set `scale` and `breaks`. It makes the function confused.")
    }

    # if color is a color mapping function
    if(is.function(color)) {
        ht_param$col = color
        if(!identical(breaks, NA)) {
            warning_wrap("`breaks` is ignored when `color` is set as a color mapping function.")
        }
    } else {
        if(identical(breaks, NA)) {
            n_col = length(color)
            if(identical(scale, "row") || identical(scale, "column")) {
                lim = max(abs(mat), na.rm = TRUE)
                ht_param$col = colorRamp2(seq(-lim, lim, length = n_col), color)
            } else {
                ht_param$col = colorRamp2(seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length = n_col), color)
            }
        } else  {
            if(length(breaks) == length(color) + 1) {
                ht_param$col = local({
                    breaks = breaks
                    color = color
                    fun = function(x) {
                        n = length(color)
                        df = data.frame(start = c(-Inf, breaks[seq_len(n)], breaks[n+1]), 
                                        end = c(breaks[1], breaks[1+seq_len(n)], Inf))
                        # tell which interval x is in
                        ind = numeric(length(x))
                        for(i in seq_along(x)) {
                            ind[i] = which(df$start <= x[i] & df$end > x[i])
                        }
                        ind = ind - 1
                        ind[ind < 1] = 1
                        ind[ind > n] = n
                        color[ind]
                    }
                    attr(fun, "breaks") = breaks
                    fun
                })
            } else if(length(breaks) == length(color)) {
                ht_param$col = colorRamp2(breaks, color)
            } else {
                n_col = length(color)
                ht_param$col  = colorRamp2(seq(min(breaks), max(breaks), length = n_col), color)
                warning_wrap("`breaks` does not have the same length as `color`. The colors are interpolated from the minimal to the maximal of `breaks`.")
            }
        }
    }
    
    if(!identical(filename, NA)) {
        warning_wrap("argument `filename` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    if(!identical(width, NA)) {
        warning_wrap("argument `width` is not supported in pheatmap -> Heatmap translation, skip it.")
    }
    
    if(!identical(height, NA)) {
        warning_wrap("argument `height` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    if(!identical(silent, FALSE)) {
        warning_wrap("argument `silent` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    ht_param$rect_gp = gpar(col = border_color)

    if(nrow(mat) > 1000 || ncol(mat) > 1000) {
        if(!is.na(border_color)) {
            warning_wrap("border color is set for the matrix with large numbers of rows or columns. You might only be able to see the border colors in the plot. Set `border_color = NA` to get rid of it.")
        }
    }

    if(!identical(cellwidth, NA)) {
        ht_param$width = ncol(mat)*unit(cellwidth, "pt")
    }

    if(!identical(cellheight, NA)) {
        ht_param$height = nrow(mat)*unit(cellheight, "pt")
    }

    if(identical(clustering_distance_rows, "correlation")) clustering_distance_rows = "pearson"
    if(identical(clustering_distance_cols, "correlation")) clustering_distance_cols = "pearson"

    ht_param$cluster_rows = cluster_rows
    ht_param$cluster_columns = cluster_cols 
    ht_param$clustering_distance_rows = clustering_distance_rows
    ht_param$clustering_distance_columns = clustering_distance_cols 
    ht_param$clustering_method_rows = clustering_method
    ht_param$clustering_method_columns = clustering_method

    if(!is.na(cutree_rows)) {
        if(inherits(cluster_rows, c("logical", "hclust", "dendrogram"))) {
            ht_param$row_split = cutree_rows
            ht_param$row_gap = unit(4, "bigpts")
            ht_param["row_title"] = list(NULL)
        }
    }
    if(!is.na(cutree_cols)) {
        if(inherits(cluster_cols, c("logical", "hclust", "dendrogram"))) {
            ht_param$column_split = cutree_cols
            ht_param$column_gap = unit(4, "bigpts")
            ht_param["column_title"] = list(NULL)
        }
    }
    
    ht_param$row_dend_width = unit(treeheight_row, "pt")
    ht_param$column_dend_height = unit(treeheight_col, "pt")

    ht_param$show_heatmap_legend = legend

    if(identical(scale, "row") || identical(scale, "column")) {
        if(identical(legend_breaks, NA)) {
            lim = quantile(abs(mat), 0.975)

            le = pretty(c(-lim, lim), n = 3)
            if(length(le) == 7 && le[1] == -3) {
                le = c(-3, -1.5, 0, 1.5, 3)
            } else if(! 0 %in% le) {
                le = c(le[1], le[1]/2, 0, le[length(le)]/2, le[length(le)])
            }
            legend_breaks = le
        }
    }
    if(!identical(legend_breaks, NA)) {
        heatmap_legend_param$at = legend_breaks
    }
    if(!identical(legend_labels, NA)) {
        heatmap_legend_param$labels = legend_labels
    }
    ht_param$heatmap_legend_param = heatmap_legend_param

    if(identical(annotation_colors, NA)) {
        annotation_colors = list()
    }
    if(!identical(annotation_col, NA)) {
        acn = rownames(annotation_col)
        mcn = colnames(mat)
        if(!is.null(acn)) {
            if(acn[1] %in% mcn) {
                if(length(union(acn, mcn)) == length(mcn)) {
                    if(!identical(acn, mcn)) {
                        warning_wrap("Column annotation has different order from matrix columns. Adjust the column annotation based on column names of the matrix.")
                    }
                    annotation_col = annotation_col[mcn, , drop = FALSE]
                }
            }
        }
        for(nm in colnames(annotation_col)) {
            if(nm %in% names(annotation_colors)) {
                if(is.null(names(annotation_colors[[nm]])) && is.numeric(annotation_col[, nm])) {
                    foo_x = annotation_col[, nm]
                    foo_n_col = length(annotation_colors[[nm]])
                    annotation_colors[[nm]] = colorRamp2(seq(min(foo_x), max(foo_x), length = foo_n_col), annotation_colors[[nm]])
                }
            }
        }
        ht_param$top_annotation = HeatmapAnnotation(df = annotation_col[, ncol(annotation_col):1, drop = FALSE], 
            col = annotation_colors, show_legend = annotation_legend,
            show_annotation_name = annotation_names_col, gp = gpar(col = border_color),
            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
            simple_anno_size = unit(10, "bigpts"), gap = unit(2, "bigpts"))
    }
    if(!identical(annotation_row, NA)) {
        arn = rownames(annotation_row)
        mrn = rownames(mat)
        if(!is.null(arn)) {
            if(arn[1] %in% mrn) {
                if(length(union(arn, mrn)) == length(mrn)) {
                    if(!identical(arn, mrn)) {
                        warning_wrap("Row annotation has different order from matrix rows. Adjust the row annotation based on row names of the matrix.")
                    }
                    annotation_row = annotation_row[mrn, , drop = FALSE]
                }
            }
        }
        for(nm in colnames(annotation_row)) {
            if(nm %in% names(annotation_colors)) {
                if(is.null(names(annotation_colors[[nm]])) && is.numeric(annotation_row[, nm])) {
                    foo_x = annotation_row[, nm]
                    foo_n_col = length(annotation_colors[[nm]])
                    annotation_colors[[nm]] = colorRamp2(seq(min(foo_x), max(foo_x), length = foo_n_col), annotation_colors[[nm]])
                }
            }
        }
        ht_param$left_annotation = rowAnnotation(df = annotation_row[, ncol(annotation_row):1, drop = FALSE], 
            col = annotation_colors, show_legend = annotation_legend,
            show_annotation_name = annotation_names_row, gp = gpar(col = border_color),
            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
            simple_anno_size = unit(10, "bigpts"), gap = unit(2, "bigpts"))
    }

    if(!identical(annotation, NA)) {
        warning_wrap("argument `annotation` is not supported in pheatmap -> Heatmap translation, skip it.")
    }

    if(identical(drop_levels, FALSE)) {
        warning_wrap("argument `drop_levels` is enfored to be TRUE, skip it.")
    }

    ht_param$show_row_names = show_rownames
    ht_param$show_column_names = show_colnames
    
    ht_param$row_names_gp = gpar(fontsize = fontsize_row)
    ht_param$column_names_gp = gpar(fontsize = fontsize_col)

    angle_col = match.arg(angle_col)[1]
    angle_col = switch(angle_col, 
                        "0" = 0,
                        "45" = 45,
                        "90" = 90,
                        "270" = 90,
                        "315" = -45)
    ht_param$column_names_rot = angle_col
    if(angle_col == 0) {
        ht_param$column_names_centered = TRUE
    }

    if(is.logical(display_numbers)) {
        if(display_numbers) {
            ht_param$layer_fun = local({
                number_format = number_format
                number_color = number_color
                fontsize_number = fontsize_number
                mat = mat
                function(j, i, x, y, w, h, fill) {
                    grid.text(sprintf(number_format, pindex(mat, i, j)), x = x, y = y, gp = gpar(col = number_color, fontsize = fontsize_number))
                }
            })
        }
    } else if(is.matrix(display_numbers)) {
        if(!identical(dim(display_numbers), dim(mat))) {
            stop_wrap("dimension of `display_numbers` should be the same as the input matrix.")
        }
        ht_param$layer_fun = local({
            number_color = number_color
            fontsize_number = fontsize_number
            mat = display_numbers
            function(j, i, x, y, w, h, fill) {
                grid.text(pindex(mat, i, j), x = x, y = y, gp = gpar(col = number_color, fontsize = fontsize_number))
            }
        })
    }
    
    if(!is.null(labels_row)) {
        ht_param$row_labels = labels_row
    }

    if(!is.null(labels_col)) {
        ht_param$column_labels = labels_col
    }

    if(!is.null(gaps_row)) {
        if(inherits(cluster_rows, c("hclust", "dendrogram"))) {
            stop_wrap("`gaps_row` should not be set when `cluster_rows` is set as a clustering object.")
        }
        if(identical(cluster_rows, TRUE)) {
            stop_wrap("`gaps_row` should not be set when `cluster_rows` is set to TRUE.")
        }
        slices = diff(c(0, gaps_row, nrow(mat)))
        ht_param$row_split = rep(seq_along(slices), times = slices)
        ht_param$row_gap = unit(4, "bigpts")
        ht_param["row_title"] = list(NULL)
    }
    if(!is.null(gaps_col)) {
        if(inherits(cluster_cols, c("hclust", "dendrogram"))) {
            stop_wrap("`gaps_col` should not be set when `cluster_cols` is set as a clustering object.")
        }
        if(identical(cluster_cols, TRUE)) {
            stop_wrap("`gaps_col` should not be set when `cluster_cols` is set to TRUE.")
        }
        slices = diff(c(0, gaps_col, ncol(mat)))
        ht_param$column_split = rep(seq_along(slices), times = slices)
        ht_param$column_gap = unit(4, "bigpts")
        ht_param["column_title"] = list(NULL)
    }

    if(!identical(clustering_callback, NA)) {
        if(!identical(ht_param$cluster_rows, FALSE)) {
            row_hclust = hclust(get_dist(mat, ht_param$clustering_distance_rows), ht_param$clustering_method_rows)
            row_hclust = clustering_callback(row_hclust, ...)
            ht_param$cluster_rows = row_hclust
        }
        if(!identical(ht_param$cluster_columns, FALSE)) {
            column_hclust = hclust(get_dist(t(mat), ht_param$clustering_distance_columns), ht_param$clustering_method_columns)
            column_hclust = clustering_callback(column_hclust, ...)
            ht_param$cluster_columns = column_hclust
        }
    }

    ### default from pheatmap
    ht_param$name = name
    ht_param$row_dend_reorder = FALSE
    ht_param$column_dend_reorder = FALSE

    if(!identical(main, NA)) {
        ht_param$column_title = main
        ht_param$column_title_gp = gpar(fontface = "bold", fontsize = 1.3*fontsize)
    }
    ht_param = c(ht_param, list(...))
    ht = do.call(Heatmap, ht_param)
    attr(ht, "translate_from") = "pheatmap"

    if(run_draw) {
        draw(ht)
    } else {
        ht
    }
}

# == title
# Compare heatmaps between pheatmap::pheatmap() and ComplexHeatmap::pheatmap()
#
# == param
# -... The same set of arguments passed to ``pheatmap::pheatmap`` and ``ComplexHeatmap::pheatmap``.
#
# == details
# The function plots two heatmaps, one by ``pheatmap::pheatmap`` and one by ``ComplexHeatmap::pheatmap``.
# Users can see the difference between the two implementations.
#
# == example
# mat = matrix(rnorm(100), 10)
# compare_pheatmap(mat)
compare_pheatmap = function(...) {
    if(!requireNamespace("pheatmap")) {
        stop_wrap("pheatmap package should be installed.")
    }
    p1 = pheatmap::pheatmap(..., silent = TRUE)$gtable
    p2 = grid.grabExpr(draw(pheatmap(...)))
    grid.newpage()
    pushViewport(viewport(x = 0, width = 0.5, y = 0, height = unit(1, "npc") - unit(1, "cm"), just = c("left", "bottom")))
    grid.draw(p1)
    popViewport()
    pushViewport(viewport(x = 0, width = 0.5, y = 1, height = unit(1, "cm"), just = c("left", "top")))
    grid.text("pheatmap::pheatmap()")
    popViewport()
    pushViewport(viewport(x = 0.5, width = 0.5, y = 0, height = unit(1, "npc") - unit(1, "cm"), just = c("left", "bottom")))
    grid.draw(p2)
    popViewport()
    pushViewport(viewport(x = 0.5, width = 0.5, y = 1, height = unit(1, "cm"), just = c("left", "top")))
    grid.text("ComplexHeatmap::pheatmap()")
    popViewport()
}

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
