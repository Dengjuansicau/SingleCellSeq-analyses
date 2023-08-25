library(velocyto.R)
library(Seurat)
epi_g <- readRDS("goat_epi.rds")
cellID <- Cells(epi_g)
cell_embeddings <- Embeddings(epi_g, reduction = "umap")
cell_embeddings_tsne <- Embeddings(epi_g, reduction = "umap")
cluster <- epi_g@meta.data$seurat_clusters
celltypes <- epi_g@meta.data$celltypes
adata<-read.loom.matrices("merged.loom")
emat <- adata$spliced
nmat <- adata$unspliced
colnames(emat) <- gsub("x","-1",colnames(emat))
colnames(emat) <-  gsub(":","_",colnames(emat))
colnames(nmat) <- gsub("x","-1",colnames(nmat))
colnames(nmat) <-  gsub(":","_",colnames(nmat))
cell.dist <- as.dist(1-armaCor(t(cell_embeddings)))
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=8)