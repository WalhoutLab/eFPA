library(pheatmap)
library(RColorBrewer)
mat = read.csv('output/heatmapTbl.csv')
boundary = read.csv('output/heatmapTbl_boundaries.csv')

heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row

#heatTbl = heatTbl / apply(heatTbl, 1,max)
# Pairwise correlation between rows (genes)
# rows.cor <- cor(t(heatTbl), use = "pairwise.complete.obs", method = "pearson")
library(stringr)
colors = read.csv('./../3_YeastJoshua_dataset/figures/simpsons_color.csv',header = F)

dev.off()
pdf('figures/boundary_heatmap.pdf',width = 24,height = 60)
labels_row0 = paste(rownames(heatTbl))
pheatmap(heatTbl, breaks = seq(-0.3,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                   "RdYlBu")))(length(seq(-0.3,1,0.001))),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         cellwidth = 12, cellheight = 5,
         cluster_rows = TRUE,cluster_cols = FALSE)

dev.off()


##### the enrichment matrix
mat = read.csv('output/heatmapTbl_pEnri.csv')
boundary = read.csv('output/heatmapTbl_boundaries_pEnri.csv')

heatTbl = -log10(mat[,2:ncol(mat)])
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row

#heatTbl = heatTbl / apply(heatTbl, 1,max)
# Pairwise correlation between rows (genes)
# rows.cor <- cor(t(heatTbl), use = "pairwise.complete.obs", method = "pearson")
library(stringr)
colors = read.csv('./../3_YeastJoshua_dataset/figures/simpsons_color.csv',header = F)

dev.off()
pdf('figures/boundary_heatmap_enrich.pdf',width = 24,height = 20)
labels_row0 = paste(rownames(heatTbl))
pheatmap(heatTbl, breaks = seq(0,15,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                       "RdYlBu")))(length(seq(0,15,0.001))),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         cellwidth = 12, cellheight = 3,
         cluster_rows = FALSE,cluster_cols = FALSE)

dev.off()

