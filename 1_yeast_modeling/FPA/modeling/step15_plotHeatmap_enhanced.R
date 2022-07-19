library(pheatmap)
library(RColorBrewer)
# the annotation of reactions
# the manual connected pathway was used; pathways with reactions less than 3 were grouped to others; isolated transporters are labeled out (inner and exchange)
mat = read.csv('output/relCorr_heatmapTbl_realDist_biomass_met.csv')
boundary = read.csv('output/heatmapTbl_boundaries_realDist_biomass_met.csv')

heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row
dev.off()
pdf('figures/boundary_heatmap_realDist_biomass_met.pdf',width = 30,height = 17.5)
# we do row-wise left nearest impute: the prediction should be same to the left-most one if no new reaction can be seen as the distance boundary increase
for (i in 1:nrow(heatTbl)){
  for (j in 1:(ncol(heatTbl))){
    if (is.na(heatTbl[i,j])){
      p1 = heatTbl[i,1:j-1]
      p1 = p1[!is.na(p1)]
      heatTbl[i,j] = p1[length(p1)]#mean(c(p1[length(p1)],p2[1]))
    }
  }
}

# reorder according to annotation table 
library(xlsx)

labels_row0 = rownames(heatTbl)
# also add the PCC values (max PCC for each row) in the annotation
library(matrixStats)
rawPCCmat = read.csv('output/PCC_titration_producingPotential.csv')
rawPCCmat = rawPCCmat[match(rownames(heatTbl),rawPCCmat$Row),]
maxPCC = rowMaxs(as.matrix(rawPCCmat[,2:ncol(rawPCCmat)]))
labels_row0 = paste(labels_row0,round(maxPCC,2))

pheatmap(heatTbl[,1:32], breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                      "RdYlBu")))(1000),
         gaps_col = 1,
         border_color = 'black',
         fontsize = 8,fontsize_row = 6, fontsize_col = 8,
         cellwidth = 12, cellheight = 6,
         cluster_rows = FALSE,cluster_cols = FALSE
         ,labels_row = labels_row0
         #,display_numbers = sigTbl[,1:23], number_color = 'Black',fontsize_number = 6
         # 05202022 we decided to not show the sig. marker since it decreases the interpretability of the figure
)

dev.off()
