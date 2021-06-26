library(pheatmap)
library(RColorBrewer)
mat = read.csv('output/subsys_annotation_for_enriched_rxns_heatmap.csv')
boundary = read.csv('output/subsys_annotation_for_enriched_rxns_colnames.csv')

heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row

dev.off()
pdf('figures/subsys_annotation_of_enriched_rxns.pdf',width = 24,height = 14)
pheatmap(heatTbl, breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                   "RdYlBu")))(1000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         cellwidth = 15, cellheight = 5,
         cluster_rows = TRUE,cluster_cols = TRUE,
         angle_col = 45)

dev.off()

