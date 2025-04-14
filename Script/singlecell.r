# jhyu
# 2025-02-27
# single_cell_RNA-seq

suppressWarnings({
  library(SeuratWrappers) 
  library(DoubletFinder)
  library(enrichplot)
  library(batchelor) 
  library(GO.db)
  library(GSEABase)
  library(org.Hs.eg.db)
  library(DOSE) 
  library(clusterProfiler)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(patchwork) 
  library(pheatmap)
  library(magrittr) 
  library(tidyverse)
  library(Seurat)
  library(dplyr)
  library(ggstatsplot)
  library(tidyr)
  library(Matrix)
  library(infercnv)
  library(tibble)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(dittoSeq)
  library(plyranges)
  library(monocle)
  library(ggsci)
  library(igraph)
  library(CellChat)
  library(NMF)
  library(ggalluvial)
})
source("./qc.r")

#1 单细胞数据质控
scRNA1 = readRDS('data_ob_v3.rds')
dir.create("1.SingleCell_QC")
#1.1 整理数据
Idents(scRNA1) <- 'orig.ident' 
scRNA1[["percent.mito"]] = PercentageFeatureSet(scRNA1,pattern = "^MT-") 
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
HB_m <- match(HB.genes, rownames(scRNA1@assays$RNA)) 
HB.genes <- rownames(scRNA1@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA1[["percent.HB"]]<-PercentageFeatureSet(scRNA1, features=HB.genes)
beforeQC_vlnplot = VlnPlot(scRNA1, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mito",
                                        "percent.HB"), 
                           ncol = 4, 
                           pt.size = 0)
ggsave("1.SingleCell_QC/BeforeQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.pdf", plot = beforeQC_vlnplot)
ggsave("1.SingleCell_QC/BeforeQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.png", plot = beforeQC_vlnplot)

#1.2 数据过滤
filter_params = c("nFeature_RNA","nCount_RNA","percent.mito")
lower_threshold = sapply( unlist(strsplit(c("NULL","NULL"," -Inf"), ",") ),
  function(x) ifelse(x == "NULL", NA, as.numeric(x)) )
if ( length(lower_threshold) != length(filter_params) ){
  stop("The lower threshold setting is not consistant with the parameters in --filters!")}
names(lower_threshold) = filter_params
upper_threshold = sapply( unlist(strsplit(c("NULL","NULL",1), ",") ),
  function(x) ifelse(x == "NULL", NA, as.numeric(x)) )
if ( length(upper_threshold) != length(filter_params) ){
    stop("The upper threshold setting is not consistant with the parameters in --filters!")}
names(upper_threshold) = filter_params
bounds_list = list()
for ( x in filter_params)  bounds_list[[x]] = c(min = unname(lower_threshold[x]), max = unname(upper_threshold[x]) )
outliers = FindOutliers(scRNA1, vars = filter_params,
                                            var.limit = bounds_list, batch = "sampleid",
                                            type = "both", cut.1 = "mean", cut.2 = "sd",
                                            n = 2, log = FALSE )
outliercells = do.call(cbind, outliers)
metric_outlier = apply(outliercells, 1, function(x) any(x == T))
scRNA1 = AddMetaData(scRNA1, metadata = metric_outlier,
                                            col.name = "is_metric_outlier")
outlier_variables = "is_metric_outlier"
is_valid_cell = !apply(FetchData(scRNA1, vars = outlier_variables), 1, function(x) any(x == T))
scRNA1 = AddMetaData(scRNA1, metadata = is_valid_cell, col.name = "is_valid")
scRNA1 = subset(scRNA1, subset = is_valid == TRUE )
counts<-GetAssayData(scRNA1,'counts')
gmt_list = unlist(strsplit("/data/database/cellranger-refdata/refdata-gex-GRCh38-2024-A/MT_genelist.gmt",",",perl = T))
gset_list <- lapply(gmt_list, function(gmtfile){
                            gset_list <- GSEABase::geneIds(GSEABase::getGmt(con=gmtfile))  
                            return(gset_list)
        })
char_vector <- unlist(gset_list[[1]])
mt <- counts[char_vector,]
mt_umi = colSums(mt)
median_value <- summary(mt_umi)["Median"]
fust = as.data.frame(mt_umi)
fust$bak = fust$mt_umi
filter_cell = rownames(fust[(fust$mt_umi < median_value*4),])
scRNA1 = scRNA1[,filter_cell]
scRNA1@meta.data$log10GenesPerUMI <- log10(scRNA1@meta.data$nFeature_RNA)/log10(scRNA1@meta.data$nCount_RNA)
obj = SplitObject(scRNA1, split.by = "orig.ident")
obj_rm=list() 
doublets_plot = list() 
for( i in names(obj)){
    print(i)
    obj[[i]] <- NormalizeData(obj[[i]])
    obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
    obj[[i]] <- ScaleData(obj[[i]])
    obj[[i]] <- RunPCA(obj[[i]])
    obj[[i]] <- RunUMAP(obj[[i]], dims = 1:30)
    obj[[i]] <- FindNeighbors(obj[[i]], dims = 1:30) %>% FindClusters(resolution = 0.3)
    tmp <- RemoveDoublets(obj[[i]], doublet.rate=0.008,pc.num=1:30)
    obj_rm[[i]] <- tmp$obj 
    doublets_plot[[i]] <- tmp$plot
  }
scRNA1 <- obj_rm[[1]] 
if(length(obj_rm) > 1) { 
  for(i in 2:length(obj_rm)) {
    scRNA1 <- merge(scRNA1, y = obj_rm[[i]]) 
  }
}
Idents(scRNA1) <- 'orig.ident' 
afterQC_vlnplot = VlnPlot(scRNA1, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mito",
                                        "percent.HB"), 
                           ncol = 4, 
                           pt.size = 0) 

ggsave("1.SingleCell_QC/afterQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.pdf",plot = afterQC_vlnplot,width = 8,height = 6)
ggsave("1.SingleCell_QC/afterQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.png", plot = afterQC_vlnplot,width = 8,height = 6)

#1.3 数据归一化与标准化
scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1,  selection.method = "vst") 
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))

#2 批次矫正和降维聚类
dir.create("2.Clustering")
scRNAlist <- SplitObject(scRNA1, split.by = "orig.ident")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))
scRNA_mnn <- RunFastMNN(object.list = scRNAlist) 
scRNA_mnn <- FindVariableFeatures(scRNA_mnn)
scRNA_mnn <- RunTSNE(scRNA_mnn, reduction = "mnn", dims = 1:12)
scRNA_mnn <- FindNeighbors(scRNA_mnn, reduction = "mnn", dims = 1:12)
scRNA_mnn <- FindClusters(scRNA_mnn,resolution = 0.8)
scRNA_mnn@meta.data$clusters = as.numeric(scRNA_mnn@meta.data$seurat_clusters)
scRNA_mnn@meta.data$clusters = as.factor(scRNA_mnn@meta.data$clusters)
Idents(scRNA_mnn) = "clusters"
p3 <- DimPlot(scRNA_mnn,  pt.size=0.1,label = T) 
p4 <- DimPlot(scRNA_mnn, split.by="orig.ident", pt.size=0.1)
p = p3 + p4 + plot_layout(guides='collect')
ggsave('2.Clustering/MNN.pdf', p, width=8, height=4)
ggsave('2.Clustering/MNN.png', p, width=8, height=4)
saveRDS(scRNA_mnn,"MNN.rds")

#3 相关性分析
dir.create("3.Correlation")
groupby = "clusters"
groupby_data = vector()
scRNA_mnn@meta.data = droplevels(scRNA_mnn@meta.data)
scRNA_mnn@meta.data[,groupby] = as.character(scRNA_mnn@meta.data[,groupby])
for (i in names(table(scRNA_mnn[["clusters"]]))) {
    sub_ob = scRNA_mnn[, rownames(scRNA_mnn@meta.data[which(scRNA_mnn@meta.data[,groupby] %in% i),])]
    normalized_data = as.matrix(sub_ob[["RNA"]]@data)
    meta.data = sub_ob@meta.data %>% tibble::rownames_to_column(var = "id")
    groupby_data = cbind(groupby_data,rowMeans(normalized_data))
}
colnames(groupby_data) = names(table(scRNA_mnn[["clusters"]]))
data = tibble::rownames_to_column(as.data.frame(groupby_data),var="GeneID")
write.table(data, file.path("3.Correlation",paste0("normalized_data_groupby_",groupby,".xls")),quote = F, row.names = F, sep = "\t")

colnames(groupby_data) = gsub('^',paste0("clusters","_"),colnames(groupby_data))
matrix<-cor(groupby_data,method="pearson")
wid<-5+1.5*log2(length(colnames(data)))
hig<-5+1.5*log2(length(colnames(data)))
coefficient = pheatmap::pheatmap(matrix,
                      display_numbers = F,
                      border_color = "white",
                      scale = "none",
                      fontsize_number=(10.0+0.0001*log2(length(colnames(data)))),
                      number_format = "%.1f",
                      fontsize_row = (10.0+0.0001*log2(length(colnames(data)))),
                      fontsize_col = (10.0+0.0001*log2(length(colnames(data)))),
                      number_color="black",
                      angle_col=45)
ggsave('3.Correlation/coefficient_heatmap.pdf', coefficient, width=8, height=4)
ggsave('3.Correlation/coefficient_heatmap.png', coefficient, width=8, height=4)

#4 marker基因鉴定
dir.create("4.Marker")
all.markers = FindAllMarkers(scRNA_mnn,min.pct = 0.25, 
                            logfc.threshold = 0, 
                            only.pos = TRUE)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(all.markers, 
            "4.Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
write.table(top10, 
            "4.Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
scRNA_mnn <- ScaleData(scRNA_mnn, features = row.names(scRNA_mnn))
heatmap_plot = DoHeatmap(object = scRNA_mnn, 
                         features = as.vector(top10$gene), 
                         group.by = "clusters", 
                         group.bar = T, 
                         size = 3) +
  theme(axis.text.y = element_text(size = 4))
ggsave("4.Marker/top10_marker_of_each_cluster_heatmap.pdf", width = 12, height = 12,
       plot = heatmap_plot)
ggsave("4.Marker/top10_marker_of_each_cluster_heatmap.png", width = 12, height = 12,
       plot = heatmap_plot)

#5 差异基因鉴定
dir.create("5.Diffexp")
contrasts = unlist(strsplit(c("group:N:T"), ":", perl = T))
numerator_subset = unlist(strsplit(contrasts[2], ",", perl = T))
denominator_subset = unlist(strsplit(contrasts[3], ",", perl = T))
cellmeta = Seurat::FetchData(scRNA_mnn, vars = contrasts[1]) %>% tibble::rownames_to_column(var = "barcode")
numerator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
    numerator_subset) %>% dplyr::pull(barcode)
denominator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
    denominator_subset) %>% dplyr::pull(barcode)
Diff_exp = FindMarkers(scRNA_mnn, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, 
                               only.pos = TRUE,
                               ident.1 = numerator_subset,
                               ident.2 = denominator_subset, 
                               group.by = contrasts[1],
                               test.use = "wilcox")
Diff_exp = Diff_exp[abs(Diff_exp$avg_log2FC) > 0.25,]
Diff_exp = Diff_exp %>% tibble::rownames_to_column(var = "gene") %>%
        dplyr::rename(pvalue = p_val, padj = p_val_adj)
numerator_means = Matrix::rowMeans(SeuratObject::GetAssayData(scRNA_mnn,
        slot = "data")[Diff_exp$gene, numerator])
denominator_means = Matrix::rowMeans(SeuratObject::GetAssayData(scRNA_mnn,
        slot = "data")[Diff_exp$gene, denominator])
Diff_exp1 = Diff_exp %>% dplyr::mutate(FoldChange = 2^avg_log2FC, baseMean = 1/2 *
        (log2(numerator_means) + log2(denominator_means))) %>%
        dplyr::rename(log2FoldChange = avg_log2FC) %>% dplyr::select(gene,
        everything()) %>% dplyr::select(-baseMean)
colnames(Diff_exp1) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange")
res_Significant = dplyr::filter(Diff_exp1, `p-value` < 0.05,
            abs(log2FoldChange) > log2(1.5))
res_Significant[which(res_Significant$log2FoldChange > 0),
        "Regulation"] <- "Up"
res_Significant[which(res_Significant$log2FoldChange < 0),
        "Regulation"] <- "Down"
colnames(res_Significant) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange","Regulation")
write.table(Diff_exp, 
            "5.Diffexp/all_diff.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(res_Significant, 
            "5.Diffexp/diff_p<0.05_FC>1.5.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

#6 富集分析
#6.1 GO富集分析
dir.create("6.enrichment")
genes_symbol <- as.character(res_Significant$gene)
eg = bitr(genes_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
ego <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
GO_dot = dotplot(ego_all,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free") 
GO_bar = barplot(ego_all,split = "ONTOLOGY")+ facet_grid(ONTOLOGY~.,scales = "free")
res_plot <- CombinePlots(list(GO_dot,GO_bar), nrow=1)
ggsave("6.enrichment/GO_results_all.pdf", plot=res_plot, width = 12,height = 10)
ggsave("6.enrichment/GO_results_all.png", plot=res_plot, width = 12,height = 10)

#6.2 GSEA富集分析
geneList= Diff_exp$avg_log2FC 
names(geneList)= toupper((Diff_exp$gene))
geneList=sort(geneList,decreasing = T)
gmtfile ='c5.go.bp.v2023.1.Hs.symbols.gmt'
geneset <- read.gmt( gmtfile )  
egmt <- GSEA(geneList, TERM2GENE=geneset, 
              minGSSize = 1,
              pvalueCutoff = 1,
              verbose=FALSE)
gsea_results_df <- egmt@result
rownames(gsea_results_df)
write.table(gsea_results_df,file = '6.enrichment/gsea_results_df.xls',quote = F,row.names = F,sep = '\t')
p = gseaplot2(egmt,geneSetID = 'KEGG_Glioma(hsa05214)',pvalue_table=T)
ggsave('6.enrichment/gsea_results_df.pdf',p,width = 8,height = 4)
ggsave('6.enrichment/gsea_results_df.png',p,width = 8,height = 4)

#6.3 KEGG富集分析
# https://davidbioinformatics.nih.gov/conversion.jsp
# http://bioinfo.org/kobas/genelist/

#7 AddModuleScore
dir.create("7.AddModuleScore")
extra_gene = read.delim("geneset.addmodeulescore.xls",check.names = FALSE)
score_outdir = "7.AddModuleScore"
if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1]="extra"}
formated_extra_gene = as.data.frame(tidyr::gather(extra_gene,key = "cluster",value = "GENE"))
match_results = CaseMatch(search = as.vector(formated_extra_gene$GENE),match = rownames(scRNA_mnn))
match_results <- Filter(function(x) length(x) > 0, match_results)
filtered_gene = formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match_results )& formated_extra_gene$GENE != ""]
if(length(filtered_gene) != 0){
filtered_gene = as.data.frame(filtered_gene)
colnames(filtered_gene) = "Gene"
write.table(filtered_gene,file.path(score_outdir,"genes_not_matched.xls"),quote = F,row.names=F)
print("有部分基因未匹配到，见：genes_not_matched.xls.")
}
formated_extra_gene = formated_extra_gene %>% 
            dplyr::filter(GENE %in% names(match_results)) %>%
            rowwise() %>%
            mutate(MATCHED = list(match_results[[GENE]])) %>%
            unnest(cols = MATCHED) %>%
            dplyr::rename(folder_suffix = cluster, gene = MATCHED)
topn_markers = formated_extra_gene
topn_markers2vis=list()
for ( clusterx in unique(topn_markers$folder_suffix) ){
    topn_markers2vis[[clusterx]] = subset(topn_markers,folder_suffix == clusterx)$gene
}
scRNA_mnn = AddModuleScore(scRNA_mnn,features=topn_markers2vis,name=names(topn_markers2vis),seed=1)
colnames(scRNA_mnn@meta.data)[(dim(scRNA_mnn[[]])[2]-length(topn_markers2vis)+1):dim(scRNA_mnn[[]])[2]] = names(topn_markers2vis)
matrix = scRNA_mnn@meta.data %>%
                    dplyr::rename( "Barcode" = "rawbc") %>%
                    dplyr::select( Barcode,"clusters",!!names(topn_markers2vis) ) %>%
                    dplyr::arrange(!!sym("clusters"))
write.table(matrix, quote = F,sep ="\t",row.names = F,file.path(score_outdir,paste0("addmodeulescore_plot.xls",collapse = "")))

for (i in colnames(matrix)[-(1:2)]){
p = ggstatsplot::ggbetweenstats(matrix,x = !!sym("clusters"),y = !!sym(i) ,
                                    plot.type = "boxviolin",
                                    results.subtitle =FALSE,
                                    messages = FALSE,
                                    pairwise.comparisons =FALSE, 
                                    mean.label.size = 0,
                                    centrality.plotting = FALSE,
                                    ylab = paste(i)) + 
              scale_color_manual(values= c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3")) +
              theme(axis.text.x = element_text(size=8,colour="black",angle = 30,vjust = 0.85,hjust = 0.75),
                    axis.text.y = element_text(size=8,colour="black"),
                    panel.grid.major =element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"))
      savenames = paste0(i,"_ggstatsplot_violin_boxplot")
      ggsave(file.path(score_outdir,paste0(savenames,".pdf")),plot = p,bg="white")
      ggsave(file.path(score_outdir,paste0(savenames,".png")),plot = p,bg="white")
}
lapply(names(topn_markers2vis),function(x){
    plot= Seurat::FeaturePlot(scRNA_mnn,features = x,
                            cols = c("grey","red"),split.by = NULL, reduction = "tsne", 
                            ncol = 1, pt.size = 0.5, order = T) +
                            theme( plot.title = element_text(hjust = 0.5), 
                                    legend.position = "right") + 
                            theme(aspect.ratio = 1/1)
    ggsave(file.path(score_outdir,paste0(x,"_score_featureplot.pdf")),plot=plot,width = 6,height = 5,limitsize = FALSE)
    ggsave(file.path(score_outdir,paste0(x,"_score_featureplot.png")),plot=plot,width = 6,height = 5,bg="white")
})

# 8 InferCNV
dir.create("8.InferCNV")
gtf= plyranges::read_gff('genes.gtf')
gene.chr = gtf %>% plyranges::filter(type == "gene" & gene_name %in% rownames(scRNA_mnn)) %>%
    as.data.frame() %>%
    dplyr::select(gene_name, seqnames, start, end) %>%
    dplyr::distinct(gene_name, .keep_all=T) %>% 
    dplyr::mutate( seqnames =seqnames)
count_mat = GetAssayData(scRNA_mnn, "counts")
cellanno = FetchData(scRNA_mnn, vars = 'new_celltype' ) %>% tibble::rownames_to_column(var = "cellbarcode")
tempdir = tempdir()
cnv_celltyping = file.path(tempdir, "cnv_celltype_group.xls")
write.table(cellanno, cnv_celltyping, sep = "\t", col.names = F,row.names = F, quote = F)
gene_order_f = file.path(tempdir, "gene_order_file.xls" )
write.table(gene.chr, gene_order_f, col.names =F, row.names =F, sep = "\t", quote =F )
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= count_mat,
                                    annotations_file= cnv_celltyping,
                                    delim="\t",
                                    gene_order_file= gene_order_f,
                                    ref_group_names=c('T_cells'))
output_dir='8.InferCNV'
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff= 0.1, 
                             analysis_mode= 'subclusters',
                             tumor_subcluster_pval=0.05,
                             hclust_method = 'ward.D2', 
                             out_dir= output_dir,
                             num_threads=2,
                             cluster_by_groups=T,
                             denoise=T,
                             HMM=T)
pdf( file.path("8.InferCNV/heatmap.pdf"), width = 18, height = 12 )
ComplexHeatmap::Heatmap(
                     t(as.matrix(infercnv_obj@expr.data)),
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_row_names =F,
                     show_column_names = F,
                     name="CNV level",
                     use_raster=TRUE,
                     raster_quality=4 )
dev.off()
scRNA_mnn[['CNV']] = CreateAssayObject(data = infercnv_obj@expr.data)
infercnv_level <- apply(as.data.frame(t(infercnv_obj@expr.data)), 1, function(x) {
          x[is.na(x)] <- 0
          return(sum(x))
        })
infercnv_level <- round(scales::rescale(infercnv_level / nrow(infercnv_obj@expr.data), c(1, 100)), 0)
infercnv_level <- infercnv_level[Cells(scRNA_mnn)]
scRNA_mnn@meta.data$cnv_level = infercnv_level
p1 = FeaturePlot(scRNA_mnn , features = "cnv_level" )
p2 = VlnPlot(scRNA_mnn , features = "cnv_level" , group.by = "orig.ident",pt.size = 0)
ggsave("8.InferCNV/featureplot.png",p1,height = 5,width = 5)
ggsave("8.InferCNV/vlnplot.png",p2,height = 4,width = 7)

#9 CellCycle
dir.create("9.CellCycle")
spe = "human"
ref.pairs <- readRDS(system.file("exdata", paste0(spe,"_cycle_markers.rds"), package="scran"))
a <- read.table("human.txt",header=F)
a=tibble::column_to_rownames(as.data.frame(a),"V1")
a[,1]=as.character(a[,1])
ref=list(
    G1=data.frame(first=a[ref.pairs$G1[,1],],second=a[ref.pairs$G1[,2],],stringsAsFactors=F),
    S=data.frame(first=a[ref.pairs$S[,1],],second=a[ref.pairs$S[,2],],stringsAsFactors=F),
    G2M=data.frame(first=a[ref.pairs$G2M[,1],],second=a[ref.pairs$G2M[,2],],stringsAsFactors=F)
)
if( length(intersect(rownames(scRNA_mnn),ref[[1]]$first))==0 ) stop("Gene names not matched. please check the species.")
sce_ob<- as.SingleCellExperiment(scRNA_mnn) #to sce
assignments <- scran::cyclone(sce_ob, ref)
scores= assignments$normalized.scores %>% mutate(barcodes=rownames(scRNA_mnn@meta.data), phase=assignments$phases) %>% select(barcodes,everything()) %>% tibble::column_to_rownames("barcodes")
colnames(scores) <- paste0("scrna","_",colnames(scores))
colnames(scores)[4] <- paste0("scrna","_CellCycle")
scRNA_mnn <- AddMetaData(scRNA_mnn, metadata=scores, col.name =colnames(scores) )
ggdim = DimPlot(object = scRNA_mnn,
               dims = c(1,2),
               reduction = "tsne",
               pt.size = 0.1,
               group.by = "scrna_CellCycle" ) +
               ggplot2::scale_colour_manual( values = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3"))
ggsave('9.CellCycle/featureplot.pdf', ggdim, width=8, height=4)
ggsave('9.CellCycle/featureplot.png', ggdim, width=8, height=4)
meta.data = scRNA_mnn@meta.data %>%
                mutate( cell_barcode = rownames(scRNA_mnn@meta.data)) %>%
                dplyr::select(cell_barcode,everything())
write.table(meta.data, file.path("9.CellCycle/", "cell_cycle_annotation_result.xls"),
            col.names =T, row.names =F,sep="\t",quote=F )

#10 Monocle2
dir.create("10.Monocle2")
cds <- as.CellDataSet(scRNA_mnn)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds = detectGenes(cds,min_expr = 1)
expressed_genes = row.names(subset(fData(cds), num_cells_expressed > 10))
clustering_DEGs = differentialGeneTest(cds[expressed_genes,],
                                       fullModelFormulaStr ="~clusters",cores = 4)
featureData(cds)@data[rownames(clustering_DEGs),"qval"]=clustering_DEGs$qval
ordering_genes <- row.names (subset(clustering_DEGs, qval < 0.01))
gbm_cds = setOrderingFilter(cds,ordering_genes = ordering_genes)
p = plot_ordering_genes(gbm_cds)
ggsave('10.Monocle2/ordering_genes.pdf', p, width=8, height=4)
ggsave('10.Monocle2/ordering_genes.png', p, width=8, height=4)
gbm_cds = reduceDimension(gbm_cds, max_components = 2, verbose = T,
                          check_duplicates = F, num_dim = 10)
gbm_cds = orderCells(gbm_cds, reverse = F)
p = plot_cell_trajectory(gbm_cds, color_by = "State", cell_size = 1.5, show_branch_points = T)+ scale_color_simpsons()
ggsave('10.Monocle2/cell_trajectory.pdf', p, width=8, height=4)
ggsave('10.Monocle2/cell_trajectory.png', p, width=8, height=4)
p = plot_cell_trajectory(gbm_cds, color_by = "clusters", 
                     cell_size = 1.5, show_branch_points = T)
ggsave('10.Monocle2/cell_trajectory.pdf', p, width=8, height=4)
ggsave('10.Monocle2/cell_trajectory.png', p, width=8, height=4)
p = plot_cell_trajectory(gbm_cds, color_by = "clusters", 
                     cell_size = 1.5, show_branch_points = T)+
     facet_wrap(~clusters)
ggsave('10.Monocle2/cell_trajectory_facet.pdf', p, width=8, height=4)
ggsave('10.Monocle2/cell_trajectory_facet.png', p, width=8, height=4)

p=plot_cell_trajectory(gbm_cds, color_by = "Pseudotime",
                     show_branch_points = F) +
         scale_colour_viridis_c(option = "inferno")
ggsave('10.Monocle2/cell_trajectory_pseudotime.pdf', p, width=8, height=4)
ggsave('10.Monocle2/cell_trajectory_pseudotime.png', p, width=8, height=4)
p = plot_complex_cell_trajectory(gbm_cds, color_by = "clusters", 
                                  show_branch_points = T, cell_size = 1, 
                                  cell_link_size = 0.3) 
ggsave('10.Monocle2/cell_trajectory_facet.pdf', p, width=8, height=4)
ggsave('10.Monocle2/cell_trajectory_facet.png', p, width=8, height=4)
p = plot_complex_cell_trajectory(gbm_cds, color_by = "clusters", 
                                  show_branch_points = T, 
                                  cell_size = 1, 
                                  cell_link_size = 0.3) + facet_wrap(~clusters)
ggsave('10.Monocle2/cell_trajectory_facet.pdf', p, width=8, height=4)
ggsave('10.Monocle2/cell_trajectory_facet.png', p, width=8, height=4)
genes <- as.factor(subset(gbm_cds@featureData@data, use_for_ordering == TRUE)$gene_short_name)
to_be_tested <- row.names(subset(fData(gbm_cds), gene_short_name %in% levels(genes)))
gbm_cds <- gbm_cds[to_be_tested, ]
varMetadata(gbm_cds)[,1] = rownames(varMetadata(gbm_cds))
gbm_cds@featureData@varMetadata[,1] = rownames(gbm_cds@featureData@varMetadata)
p <- plot_pseudotime_heatmap(gbm_cds, cores = 1, 
                             cluster_rows = T, num_clusters = 4, 
                             show_rownames = F, return_heatmap = T)
ggsave('10.Monocle2/pseudotime_heatmap.pdf', p$ph_res, width=8, height=4)
ggsave('10.Monocle2/pseudotime_heatmap.png', p$ph_res, width=8, height=4)
gene_clusters <- cutree(p$tree_row, k = 4)
gene_clustering <- data.frame(gene_clusters)
gene_clustering[, 1] <- as.character(gene_clustering[, 1])
colnames(gene_clustering) <- "gene_module"
BEAM_res <- BEAM(gbm_cds, branch_point = 2, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
p = plot_genes_branched_heatmap(gbm_cds[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 2,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = T)
ggsave("10.Monocle2/pseudotime_heatmap_branchtime.pdf", plot = p$ph_res,
      height =7,width = 7, bg="white")
ggsave("pseudotime_heatmap_branchtime.png", plot = p$ph_res,
      height =7,width = 7, dpi = 1000,bg="white")
p$heatmap_matrix_ori[1:5,1:10]
gene_clusters <- cutree(p$ph_res$tree_row, k = 4)
gene_clustering <- data.frame(gene_clusters)
gene_clustering[, 1] <- as.character(gene_clustering[, 1])
colnames(gene_clustering) <- "gene_module"
to_be_tested_sub <- row.names(subset(fData(gbm_cds), 
                                     gene_short_name %in% c("ISG15","PARK7", "GLUL")))
p = plot_genes_jitter(gbm_cds[to_be_tested_sub,], grouping = "State", 
                  min_expr = 0.1,color_by = "State",cell_size = 1)+ 
           scale_color_simpsons()
ggsave('10.Monocle2/genes_jitter.pdf', p, width=8, height=4)
ggsave('10.Monocle2/genes_jitter.png', p, width=8, height=4)
p = plot_genes_in_pseudotime(gbm_cds[to_be_tested_sub,],
                         color_by = "clusters", cell_size = 1, ncol = 1) +scale_color_simpsons()
ggsave('10.Monocle2/genes_in_pseudotime.pdf', p, width=8, height=4)
ggsave('10.Monocle2/genes_in_pseudotime.png', p, width=8, height=4)
branchpoint = 2
new_cds <- buildBranchCellDataSet(gbm_cds[to_be_tested_sub,],
                                  branch_point = branchpoint, 
                                  progenitor_method = "duplicate")
cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
branch_labels <- c(paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")), 
                   paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-")))
p = plot_genes_branched_pseudotime(gbm_cds[to_be_tested_sub,], 
                               color_by = "clusters", 
                               branch_point = branchpoint, 
                               cell_size = 1, ncol = 1, 
                               branch_labels = branch_labels) + scale_color_simpsons()
ggsave('10.Monocle2/genes_branched_pseudotime.pdf', p, width=8, height=4)
ggsave('10.Monocle2/genes_branched_pseudotime.png', p, width=8, height=4)
p <- plot_cell_trajectory(gbm_cds, markers = "GLUL",
                          use_color_gradient = T, show_branch_points = F, 
                          show_tree = F, cell_size = 1.5) + 
  theme(legend.text = element_text(size = 10)) + 
  scale_color_gradientn(colours = c("grey", "yellow", "red"))
ggsave('10.Monocle2/genes_branched_pseudotime.pdf', p, width=8, height=4)
ggsave('10.Monocle2/genes_branched_pseudotime.png', p, width=8, height=4)

#11 cellchat
#11.1 整体cellchat进行细胞通讯分析
dir.create("11.cellchat")
scRNA_mnn <- SetIdent(scRNA_mnn, value = "new_celltype")
cellchat <- createCellChat(object = scRNA_mnn, group.by = "new_celltype")
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
showDatabaseCategory(CellChatDB.human)
cellchat <- CellChat::subsetData(cellchat) 
cellchat@data.signaling[1:10,1:10]
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat,slot.name = "netP")
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("11.cellchat/cellchat_circle_count.pdf",width = 8,height = 4)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("11.cellchat/cellchat_circle_weight.pdf",width = 8,height = 4)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat@netP$pathways
pathways.show <- c("MIF") 
levels(cellchat@meta$celltype)
pdf("11.cellchat/cellchat_circle_aggregate.pdf",width = 8,height = 4)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
pdf("11.cellchat/cellchat_circle_aggregate_hierarchy.pdf",width = 8,height = 4)
vertex.receiver = seq(1,3)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()
pdf("11.cellchat/cellchat_circle_aggregate_chord.pdf",width = 8,height = 4)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
pdf("11.cellchat/cellchat_heatmap_count.pdf",width = 8,height = 4)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds")
dev.off()
pdf("11.cellchat/cellchat_heatmap_weight.pdf",width = 8,height = 4)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
dev.off()
pdf("11.cellchat/cellchat_heatmap_pathways.pdf",width = 8,height = 4)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("11.cellchat/cellchat_heatmap_pathways_centrality.pdf",width = 8,height = 4)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show)
dev.off()
pdf("11.cellchat/cellchat_heatmap_pathways_centrality.pdf",width = 8,height = 4)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()
pdf("11.cellchat/cellchat_heatmap_pathways_centrality.pdf",width = 8,height = 4)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.off()
pdf("11.cellchat/cellchat_bubble.pdf",width = 8,height = 4)
netVisual_bubble(cellchat, remove.isolate = FALSE, return.data = FALSE)
dev.off()
pdf("11.cellchat/cellchat_chord_gene.pdf",width = 8,height = 4)
netVisual_chord_gene(cellchat, lab.cex = 0.8,legend.pos.y = 50,legend.pos.x = 50)
dev.off()

#11.2 分别对各组cellchat进行细胞通讯分析
scRNA_mnn@meta.data %>% group_by(orig.ident,new_celltype) %>% summarise(n = n())
scRNA_mnn = subset(scRNA_mnn, 
                   new_celltype %in% c("EC1","EC3","EC4","FC1","FC2","FC3"))
scRNA_mnn@meta.data$celltype = droplevels(scRNA_mnn@meta.data$celltype)
sp <- SplitObject(scRNA_mnn, split.by = "orig.ident")
cellchat_A <- createCellChat(object = sp[["HNC01TIL"]], 
                             group.by = "new_celltype")
cellchat_B <- createCellChat(object = sp[["Tonsil2"]], 
                             group.by = "new_celltype")
cellchat_A@idents = droplevels(cellchat_A@idents)
cellchat_B@idents = droplevels(cellchat_B@idents)
cellchat_list= list()
for (i in c(cellchat_A,cellchat_B)){
  cellchat = i
  name = unique(i@meta$orig.ident)
  CellChatDB <- CellChatDB.human 
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <- CellChat::subsetData(cellchat) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat_list[[name]] = cellchat
}
cellchat <- mergeCellChat(cellchat_list , add.names = names(cellchat_list))
pdf("11.cellchat/cellchat_compare_interaction.pdf",width = 8,height = 4)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()
pdf("11.cellchat/cellchat_compare_interaction_diff.pdf",width = 8,height = 4)
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(2,1))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()
pdf("11.cellchat/cellchat_compare_interaction_heatmap.pdf",width = 8,height = 4)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1+gg2
dev.off()

# Uphyloplot2
#https://github.com/harbourlab/UPhyloplot2


# scvelo
# https://github.com/velocyto-team/velocyto.py







