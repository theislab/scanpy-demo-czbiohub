#date created: 10/10/19
#author: Maren Buettner
#place: Institute of Computational Biology, Helmholtz Centre Munich
#purpose: This is a template to read in anndata objects in R

library(rhdf5)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)

#set paths and dates (please adjust accordingly)
figure_path <- './figures/'
file_path <- './table/'
f_path <- './data/adata.h5ad'
today <- format(Sys.Date(), '%y%m%d')

#Read data
adata <- h5read(f_path,'/',compoundAsDataFrame=FALSE)

#get genes/cell ID
barcodes <- unlist(adata$obs$index)
genes <- unlist(adata$var$index)

data_raw <- adata$raw.X$data
index_raw <- adata$raw.X$indices
ptr_raw <- adata$raw.X$indptr
sparse_mat <- sparseMatrix(p = as.numeric(ptr_raw), 
                           x=as.numeric(data_raw), 
                           i = as.numeric(index_raw)+1)
genes_raw <- adata$raw.var$index

#compute n_genes, n_counts
n_genes <- colSums(sparse_mat>0)
n_counts <- colSums(sparse_mat)

#set levels
sample_levels <- c( "sample","levels", "as", "in", "uns$sample_categories")
cell_type_levels <- c("cell", "types")

#get attributes
cellData <- data.frame(
  sample=factor(unlist(adata$uns$sample_categories)[unlist(adata$obs$sample)+1], 
                levels=sample_levels),
  n_genes = n_genes,
  n_counts = n_counts,
  cell_type = factor(unlist(adata$uns$cell_type_categories)[unlist(adata$obs$cell_type+1)], 
                     levels=cell_type_levels)
) 




#create boxplot 
ggplot(cellData, aes(cell_type, n_genes, fill=sample)) +geom_boxplot() +
  scale_fill_manual(values = c('red','yellow',  'green')) + theme_bw()
ggplot(cellData, aes(cell_type, n_counts, fill=sample)) +geom_boxplot() +
  scale_fill_manual(values = c('red','yellow',  'green')) + theme_bw()


#############
# Differential expression test
#############