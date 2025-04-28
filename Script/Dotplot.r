library(ggplot2)
library(Seurat)

seurat_ob = readRDS('seurat.rds')
genes = read.delim('genelist.txt')
dot_plot = DotPlot(object = seurat_ob, features = as.vector(genes[,1]), cols=c('grey','#9d0142') ) + coord_flip() + guides(color = guide_colorbar(order = 1, title = "Average Expression"))
ggsave('dotplot.pdf', plot = dot_plot)
ggsave('dotplot.png', plot = dot_plot)
