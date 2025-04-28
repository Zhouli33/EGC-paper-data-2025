library(Seurat)
library(dplyr)
library(tibble)


seurat_ob = readRDS('seurat.rds')
groupby = 'new_celltype'
groupby_data = vector()
for (i in names(table(seurat_ob[[groupby]]))  ) {
    sub_ob = SubsetData(seurat_ob, subset.name= groupby,accept.value=i)
    normalized_data = as.matrix(sub_ob[['RNA']]@data)
    meta.data = sub_ob@meta.data %>% tibble::rownames_to_column(var = "id")
    groupby_data = cbind(groupby_data,rowMeans(normalized_data))
}
colnames(groupby_data) = names(table(seurat_ob[[groupby]]))
colnames(groupby_data) = gsub('^',paste0(groupby,"_"),colnames(groupby_data))
matrix = cor(groupby_data, method="pearson")
matrix =rownames_to_column(as.data.frame(matrix),var="person")
write.table(matrix, "Correlation.xls"),quote = F, row.names = F, sep = "\t")


