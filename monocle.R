library(monocle)
library(Seurat)
library(ggplot2)
library(dplyr)
rds <- readRDS(file="celltypes.rds")
seurat_object <- rds
data <- as(as.matrix(seurat_object@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds<- detectGenes(cds, min_expr = 0.1)
fData(cds)$use_for_ordering <- fData(cds)$num_cells_expressed > 0.0025 * ncol(cds)
cds <- reduceDimension(cds, max_components = 2, norm_method = 'log', num_dim = 20, reduction_method = 'tSNE', verbose = TRUE)
cds <- clusterCells(cds)
clustering_DEG_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~seurat_clusters', cores = 8)
my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
cds <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
GM_state <- function(cds){
if (length(unique(pData(cds)$State)) > 1){
E45_counts <- table(pData(cds)$State, pData(cds)$orig.ident)[,"E45"]
return(as.numeric(names(E45_counts)[which
(E45_counts == max(E45_counts))]))
} else {
return (1)
}
}
cds <- orderCells(cds, root_state = GM_state(cds))
save(cds,file=paste(celltypes[[i]],sep="_","cds_order.Rdata"))
pdf("celltypes.pdf")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "orig.ident")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "seurat_clusters")
dev.off()
my_pseudotime_de <- differentialGeneTest(cds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
c<-subset(my_pseudotime_de, qval < 10^-5)
write.table(c[order(c[,4]),],file="celltypes_pseudotime_deg.xls")

####plot genes###
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$gene_short_name
pdf("celltypes_genes.pdf")
plot_genes_in_pseudotime(cds[my_pseudotime_gene,],color_by = "orig.ident")
dev.off()
pdf("celltypes_heatmap.pdf")
plot_pseudotime_heatmap(cds[row.names(c),], num_clusters = 4,use_gene_short_name = TRUE,show_rownames = F,cores=5)
dev.off()
