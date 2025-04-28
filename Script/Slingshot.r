suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("slingshot"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("uwot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("grDevices"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("tradeSeq"))


sce <- readRDS('sce.rds')
groupby <- "clusters"
pointsize <- 1
reduct <- "UMAP"
new_celltype_pal = sce[['clusters_col']]
colData(sce)[[groupby]] <- factor( colData(sce)[[groupby]] )
aimCurve = c(1,2)
groupby_clust <- colData(sce)[[groupby]]
levels(groupby_clust) <- c(1:length(levels(colData(sce)[[groupby]])))
pdf('Clustering_slingshot_with_1_curve.pdf',width = 8)
par(mai=c(1,1,1,2))
plot(reducedDims(sce)[[reduct]], col = colData(sce)[[paste0( groupby,"_col" )]], pch = 16, asp = 1, cex = pointsize)
xy=par("usr")
for (i in paste0("curve",aimCurve)) {
  lines(SlingshotDataSet(sce)@curves[[i]], lwd = 2, col = 'black')
  legend(x=xy[2]+xinch(0.2), y=xy[4], xpd = TRUE,        #图例位置
	legend = names(new_celltype_pal),        #图例内容
	col = new_celltype_pal,             #图例颜色
	pch = 19 )    
}
dev.off()

pt <- slingPseudotime(sce)
gene <- read.delim('genes.txt', sep = "\t")
gene[, 1] <- factor(as.character(gene[, 1]))
topgenes <- Seurat::CaseMatch(search = levels(gene[, 1]), match = rownames(sce))
print("Start to fitGAM!")
set.seed(123)
sce0 = fitGAM(counts=as.matrix(assays(sce)$counts), sds=SlingshotDataSet(sce), pseudotime = pt, cellWeights = slingCurveWeights(sce), genes = topgenes, verbose = T )


colData(sce0)=cbind(colData(sce0),colData(sce))
gs <- lapply(topgenes, function(x) plotSmoothers(sce0[topgenes,], assays(sce)$counts[topgenes,], gene = x, lwd = 1, 
							alpha = 0.5, pointCol = groupby , border = TRUE) + labs(title = x) )
i <- 1
while (i <=  length(gs)){
  gs_sub <- gs[i]
  gs_sub <- gs_sub[which(!sapply(gs_sub, is.null))]

  pdf(paste0("pseudotime_genes_Smoothers_plot",i,".pdf", collapse = "_"),
				  width = 5, height = length(gs_sub)*2.5)
  grid.arrange(grobs = gs_sub, ncol=1) 
  dev.off()
  i=i+1
}


