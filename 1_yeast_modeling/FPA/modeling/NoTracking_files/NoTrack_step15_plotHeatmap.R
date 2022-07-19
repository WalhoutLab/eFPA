library(pheatmap)
library(RColorBrewer)
library(matrixStats)
# the annotation of reactions
# the manual connected pathway was used; pathways with reactions less than 3 were grouped to others; isolated transporters are labeled out (inner and exchange)
mat_prod = read.csv('output/PCC_titration_producingPotential.csv')
#mat_cons = read.csv('output/PCC_titration_consumingPotential.csv')
boundary = c('localFPA',as.character(seq(0,40,0.5)))
sig_prod = read.csv('output/producing_sig_rxns.csv')
#sig_cons = read.csv('output/consuming_sig_rxns.csv')

mat = mat_prod
heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary
rownames(heatTbl) = mat$Row
# filter roughly
heatTbl = heatTbl[rowMaxs(abs(as.matrix(heatTbl)))>0,]
# GTP was the only one that was left out for analysis, since it present all-same rFP (1). The reason is that 
# GTP can only be produced by reaction r_0800 (ATP--> GTP) that contains only side metabolites; therefore, 
# this reaction is excluded from the distance analysis and have inf distance to all reaction. In this case, 
# our algorithm will assign inf distance for GTP demand reaction such that nothing can be integrated. We could 
# manually allow it extend to r_0800, however, it would not go beyond this reaction because of the distance rule.
# therefore, network analysis is not meaningful anyway for GTP in our current framework. so we skip it.
dev.off()
pdf('figures/NoTrack_boundary_heatmap_metabolites_producing_all_non-zero.pdf',width = 30,height = 17.5)
labels_row0 = rownames(heatTbl)
# also add the PCC values (max PCC for each row) in the annotation
maxPCC1 = rowMaxs((as.matrix(heatTbl)))
maxPCC2 = rowMins((as.matrix(heatTbl)))
maxPCC = maxPCC1
maxPCC[abs(maxPCC2)>abs(maxPCC1)] = maxPCC2[abs(maxPCC2)>abs(maxPCC1)]
labels_row0 = paste(labels_row0,round(maxPCC,2))
heatTbl = t(scale(t(heatTbl),center = F, scale = rowMaxs(abs((as.matrix(heatTbl))))))
pheatmap(heatTbl, breaks = seq(-1,1,0.002),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                     "RdBu")))(1000),
         border_color = 'black',
         fontsize = 8,fontsize_row = 6, fontsize_col = 8,
         cellwidth = 6, cellheight = 6,
         cluster_rows = T,cluster_cols = FALSE
         ,labels_row = labels_row0
)

dev.off()

mat = mat_prod
library(stringr)
sig = str_remove(sig_prod$rxn1,'^NewMet_')
heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary
rownames(heatTbl) = mat$Row
# filter roughly
heatTbl = heatTbl[sig,]
dev.off()
pdf('figures/NoTrack_boundary_heatmap_metabolites_producing.pdf',width = 30,height = 17.5)
labels_row0 = rownames(heatTbl)
# also add the PCC values (max PCC for each row) in the annotation
maxPCC = rowMaxs(as.matrix(heatTbl))
labels_row0 = paste(labels_row0,round(maxPCC,2))
heatTbl = t(scale(t(heatTbl),center = F, scale = rowMaxs((as.matrix(heatTbl)))))
pheatmap(heatTbl, breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                    "RdYlBu")))(1000),
         border_color = 'black',
         fontsize = 8,fontsize_row = 6, fontsize_col = 8,
         cellwidth = 6, cellheight = 6,
         cluster_rows = T,cluster_cols = FALSE
         ,labels_row = labels_row0
)

dev.off()


# mat = mat_cons
# sig = sig_cons$rxn1
# heatTbl = mat[,2:ncol(mat)]
# colnames(heatTbl) = boundary
# rownames(heatTbl) = mat$Row
# # filter roughly
# heatTbl = heatTbl[rowMaxs(abs(as.matrix(heatTbl)))>0,]
# dev.off()
# pdf('figures/NoTrack_boundary_heatmap_metabolites_consuming_all_non-zero2.pdf',width = 30,height = 17.5)
# labels_row0 = rownames(heatTbl)
# # also add the PCC values (max PCC for each row) in the annotation
# maxPCC1 = rowMaxs((as.matrix(heatTbl)))
# maxPCC2 = rowMins((as.matrix(heatTbl)))
# maxPCC = maxPCC1
# maxPCC[abs(maxPCC2)>abs(maxPCC1)] = maxPCC2
# labels_row0 = paste(labels_row0,round(maxPCC,2))
# heatTbl = t(scale(t(heatTbl),center = F, scale = rowMaxs(abs(as.matrix(heatTbl)))))
# pheatmap(heatTbl, breaks = seq(-1,1,0.002),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                      "RdBu")))(1000),
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 6, fontsize_col = 8,
#          cellwidth = 6, cellheight = 6,
#          cluster_rows = T,cluster_cols = FALSE
#          ,labels_row = labels_row0
# )
# 
# dev.off()

