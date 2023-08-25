library(circlize)
library(ComplexHeatmap)
library(dplyr)
hvg <- read.table("DEGs.txt",header=F,sep="\t")
mat <- mat[hvg$V1, ]
col = colorRamp2(c(-2, 0, 2), c("#0033cc", "white", "#990000"))
pdf("top50_heatmap.pdf")
Heatmap(mat, 
        col=col,
        use_raster = T,
        cluster_rows = F,
        cluster_columns = T,
        show_column_names = F,
        show_row_names = T,
        column_split = meta,
column_title_gp = gpar(fill = type_cols, font = 1:3))
dev.off()