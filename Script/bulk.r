library(tidyverse)
suppressWarnings(library(BayesPrism))

# 1. Quality Control ####

pdf(paste0(dir.out,"plot.cor.phi.state.pdf"), height=5, width=5)
plot.cor.phi(input=sc.dat,
             input.labels=cell.state.labels,
             title="Cell State Correlation",
             cexRow=0.2, cexCol=0.2,
             margins=c(2,2))
dev.off()

pdf(paste0(dir.out,"plot.cor.phi.type.pdf"), height=5, width=5)
plot.cor.phi(input=sc.dat,
             input.labels=cell.type.labels,
             title="Cell Type Correlation",
             cexRow=0.5, cexCol=0.5)
dev.off()

pdf(paste0(dir.out,"plot.cor.phi.scRNA.outlier.pdf"), height=5, width=5)
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, # Ensure column names are gene symbols or ENSMEBL IDs
  cell.type.labels=cell.type.labels,
  species="hs", # Currently supports human(hs) and mouse(mm) annotations
  return.raw=TRUE # Returns data used for plotting
)
dev.off()

pdf(paste0(dir.out,"plot.cor.phi.bulk.outlier.pdf"), height=5, width=5)
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat, # Ensure column names are gene symbols or ENSMEBL IDs
  sc.input=sc.dat, # Ensure column names are gene symbols or ENSMEBL IDs
  cell.type.labels=cell.type.labels,
  species="hs", # Currently supports human(hs) and mouse(mm) annotations
  return.raw=TRUE
)
dev.off()

sc.dat.filtered <- cleanup.genes(input=sc.dat,
                                 input.type="count.matrix",
                                 species="hs",
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"),
                                 exp.cells=5)

pdf(paste0(dir.out,"plot.cor.phi.plot.bulk.vs.sc.pdf"), height=5, width=5)
plot.bulk.vs.sc(sc.input=sc.dat.filtered,
                bulk.input=bk.dat)
dev.off()

sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
                                       gene.type="protein_coding")

# 2. Build Objects ####

myPrism <- new.prism(reference=sc.dat.filtered.pc,
                     mixture=bk.dat,
                     input.type="count.matrix",
                     cell.type.labels=cell.type.labels,
                     cell.state.labels=cell.state.labels,
                     key="Epithelial_cells",
                     outlier.cut=0.01,
                     outlier.fraction=0.1)

bp.res <- run.prism(prism=myPrism, n.cores=40) 

# 3. Extract Tables --------------------------------------------------------------------

# Extract posterior mean of cell type fraction theta
theta.main <- get.fraction(bp=bp.res,
                           which.theta="final",
                           state.or.type="type")

theta.sub <- get.fraction(bp=bp.res,
                          which.theta="first",
                          state.or.type="state")

# 4. Calculate Proportions -------------------------------------------------------
# [Use the total number of major cell types as the denominator]

# Convert to long format
df_long <- as.data.frame(theta.sub) %>% 
  rownames_to_column("sample.ID") %>%
  pivot_longer(cols=-c("sample.ID"), values_to="value", names_to="cell.subtype")


# Add cell type column
df_long <- left_join(df_long, cell.2labels.1, by="cell.subtype")

# Calculate proportions
df_long <- df_long %>% 
  group_by(celltype, sample.ID) %>%
  mutate(total=sum(value),
         proportion=value/total) %>%
  ungroup()

df_wide <- df_long[,c("sample.ID", "cell.subtype", "proportion")] %>%
  pivot_wider(id_cols="sample.ID", names_from="cell.subtype", values_from="proportion") %>%
  column_to_rownames("sample.ID")