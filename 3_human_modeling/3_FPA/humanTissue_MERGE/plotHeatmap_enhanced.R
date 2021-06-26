library(pheatmap)
library(RColorBrewer)
################plot1: wtdDist and rel corr##########
mat = read.csv('output/relCorr_heatmapTbl_wtdDist.csv')
sig = read.csv('output/heatmapTbl_sigLabel.csv')
boundary = read.csv('output/heatmapTbl_boundaries.csv')

heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row
sigTbl = ifelse(sig[,2:ncol(mat)]==1, '+','')

#heatTbl = heatTbl / apply(heatTbl, 1,max)
# Pairwise correlation between rows (genes)
# rows.cor <- cor(t(heatTbl), use = "pairwise.complete.obs", method = "pearson")
library(stringr)
annotation = read.csv('manualPathwayAnnotation.csv')
colors = read.csv('figures/simpsons_color.csv',header = F)
#annotationList = table(annotation$pathway_major)
#annotationList = sort(annotationList,decreasing = T)
annotationList = c('Purine metabolism','Phenylalanine, tyrosine and tryptophan biosynthesis','Histidine metabolism','Lysine biosynthesis',
                   'Arginine biosynthesis','Methionine, folate and sulfur metabolism','Proline biosynthesis','Pyrimidine metabolism',
                   'central carbon metabolism','Other amino acids', 'transporter [exchange]',
                   'Fructose, mannose, sucrose, trahelose and glucan metabolism','FA and PL metabolism',
                   'transporter [inner]','Others')
colorList = colors$V1[2:(1+length(annotationList))]
#names(colorList) = c(names(annotationList)[2:(length(annotationList)-1)],names(annotationList)[1],names(annotationList)[(length(annotationList))])
names(colorList) = annotationList
colorList = list(pathway_major = colorList)
rownames(annotation) = annotation$rxn
dev.off()
pdf('figures/boundary_heatmap.pdf',width = 24,height = 14)
labels_row0 = paste(rownames(heatTbl),annotation[rownames(heatTbl),'pathway'])
pheatmap(heatTbl[,1:ncol(heatTbl)], breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                   "RdYlBu")))(1000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList, labels_row = labels_row0,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         cellwidth = 12, cellheight = 5,
         cluster_rows = TRUE,cluster_cols = FALSE
         ,display_numbers = sigTbl, number_color = 'Black',fontsize_number = 5
         )

dev.off()

################plot2: annotate the effective bound (mmIRS) ###########
mat = read.csv('output/mmIRS_wtdDist.csv')
sig = read.csv('output/heatmapTbl_sigLabel.csv')
boundary = read.csv('output/heatmapTbl_boundaries.csv')

heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row
sigTbl = ifelse(sig[,2:ncol(mat)]==1, '+','')

dev.off()
pdf('figures/boundary_heatmap_mmIRS_annotation.pdf',width = 24,height = 14)
labels_row0 = paste(rownames(heatTbl),annotation[rownames(heatTbl),'pathway'])
pheatmap(heatTbl[,1:ncol(heatTbl)], breaks = seq(0,10,1),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                   "RdYlBu")))(10),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList, labels_row = labels_row0,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         cellwidth = 12, cellheight = 5,
         cluster_rows = TRUE,cluster_cols = FALSE
         ,display_numbers = TRUE, number_color = 'Black',fontsize_number = 5, number_format = "%.1f"
)

dev.off()

################plot3: realDist and rel corr##########
mat = read.csv('output/relCorr_heatmapTbl_realDist.csv')
sig = read.csv('output/heatmapTbl_sigLabel_realDist.csv')
boundary = read.csv('output/heatmapTbl_boundaries_realDist.csv')

heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row
sigTbl = ifelse(sig[,2:ncol(mat)]==1, '+','')
dev.off()
pdf('figures/boundary_heatmap_realDist.pdf',width = 24,height = 14)
# we do row-wise left nearest impute: the prediction should be same to the left-most one if no new reaction can be seen as the distance boundary increase
for (i in 1:nrow(heatTbl)){
  for (j in 1:(ncol(heatTbl))){
    if (is.na(heatTbl[i,j])){
      p1 = heatTbl[i,1:j-1]
      p1 = p1[!is.na(p1)]
      # p2 = heatTbl[i,(j+1):ncol(heatTbl)]
      # p2 = p2[!is.na(p2)]
      s1 = sigTbl[i,1:j-1]
      s1 = s1[!is.na(p1)]
      heatTbl[i,j] = p1[length(p1)]#mean(c(p1[length(p1)],p2[1]))
      sigTbl[i,j] = s1[length(s1)]#mean(c(p1[length(p1)],p2[1]))
    }
  }
}
#heatTbl[is.na(heatTbl)] = 0
# reorder according to annotation table 
library(xlsx)
annTbl = read.xlsx('./output/predictionMechanism_annotation.xlsx','Sheet1',header = T)
rownames(sigTbl) = rownames(heatTbl)
heatTbl = heatTbl[annTbl$rxn,]
sigTbl = sigTbl[annTbl$rxn,]
labels_row0 = paste(rownames(heatTbl),annotation[rownames(heatTbl),'pathway'])
pheatmap(heatTbl[,1:38], breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                      "RdYlBu")))(1000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList, labels_row = labels_row0,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         cellwidth = 12, cellheight = 5,
         cluster_rows = FALSE,cluster_cols = FALSE
         ,display_numbers = sigTbl[,1:38], number_color = 'Black',fontsize_number = 5
)

dev.off()












########the coexpression of fluxes########
# the coexpression of fluxes
FluxMat = read.csv('output/relativeFluxTable.csv',row.names = 1)
sigInModel = mat$Row
heatTbl = abs(FluxMat)
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(heatTbl), use = "pairwise.complete.obs", method = "pearson")
# prepare the new annotation table 
colorList$isPredicted = c('#FD8D3CFF','#370335FF')
names(colorList$isPredicted) = c('Yes','No')
colorList$isCorrelated = c('#FD8D3CFF','#370335FF')
names(colorList$isCorrelated) = c('Yes','No')
annotation2 = data.frame(row.names = colnames(rows.cor), isPredicted = ifelse(colnames(rows.cor) %in% sigInModel,'Yes','No'),
                         isCorrelated = ifelse(colnames(rows.cor) %in% mat$Row[sigTbl[,boundary == 'expression only'] == '+'], 'Yes','No'))
annotation2$pathway_major = annotation[rownames(annotation2), 'pathway_major']
annotation2 = annotation2[,c('pathway_major','isCorrelated','isPredicted')]
pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                    "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,
         annotation_col = annotation2,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 2, fontsize_col = 4,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)

dev.off()
labels_col = paste(colnames(rows.cor),annotation[colnames(rows.cor),'pathway'])
pdf('figures/co-flux.pdf',width = 24,height = 14)
pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,
         annotation_col = annotation2,labels_col = labels_col,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 2, fontsize_col = 4,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off() 

########the coexpression of proteins(rxns)########
# it indicates co-expression is not always an indication of flux; sometimes even pathway-level coexpression fails, 
# here the topleft is TCA cycle that fails (need pathway annotation, the cluster is not only TCA) (but it is still weak pos cor; and it seems FPA rescued some)
# the coexpression of fluxes

# find some example that are co-expressed but not pathway-coexpressed (especially as isolated rxns), if they dont correlate with flux, it will be a 
# good argument for the pathway-level expression [[ try to check noncorrelated rxns in the same module but belongs to an unrelated pathway]]

ExpMat = read.csv('output/relativeExpressionTable.csv',row.names = 1)
# Pairwise correlation between rows (genes)
rownames(mat) = mat$Row
# load pearson reference
pearson_r = read.csv('output/flux_expression_pearson_correlation_rel2rel.csv',row.names = 1)
annotation2$pearson_r = pearson_r[rownames(annotation2),1]
colorList$pearson_r = c('#00FF00','#ff0000')
names(colorList$pearson_r) = c('-1','1')
rows.cor.exp <- cor(t(ExpMat), use = "pairwise.complete.obs", method = "pearson")
labels_col2 = paste(colnames(rows.cor.exp),annotation[colnames(rows.cor.exp),'pathway'])
pheatmap(rows.cor.exp, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,
         annotation_col = annotation2[,c("pathway_major",'pearson_r','isCorrelated','isPredicted')],labels_col = labels_col2,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 3, fontsize_col = 7,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off()
pdf('figures/co-expression.pdf',width = 24,height = 20)
pheatmap(rows.cor.exp, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                          "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,
         annotation_col = annotation2[,c("pathway_major",'pearson_r','isCorrelated','isPredicted')],labels_col = labels_col2,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 3, fontsize_col = 7,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off()
############ doesnt work -- nothing explained #####
# plus CBF1 ==> this is a glycine regulator at the least ==> correation is not consistent with target 
# CBF1_expression = read.csv('output/CBF1_expression.csv',row.names = 1)
# hist(cor(t(CBF1_expression),t(heatTbl)))
# corV = t(cor(t(CBF1_expression),t(heatTbl)))
# annotation2[rownames(corV),'CBF1'] = corV[,1]
# colorList$CBF1 = c('#00FF00','#ff0000')
# names(colorList$CBF1) = c('-0.65','0.65')
# CBF1_rxn = read.csv('output/CBF1_targetRxn.csv')
# annotation2$isCBF1target = ifelse(colnames(rows.cor) %in% CBF1_rxn$Var1, 'Yes','No')
# colorList$isCBF1target = c('#FD8D3CFF','#370335FF')
# names(colorList$isCBF1target) = c('Yes','No')
# 
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2,
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)

# plus Tup1
# Tup1_expression = read.csv('output/Tup1_expression.csv',row.names = 1)
# hist(cor(t(Tup1_expression),t(heatTbl)))
# corV = t(cor(t(Tup1_expression),t(heatTbl)))
# annotation2[rownames(corV),'Tup1'] = corV[,1]
# colorList$Tup1 = c('#00FF00','#ff0000')
# names(colorList$Tup1) = c('-0.65','0.65')
# Tup1_rxn = read.csv('output/Tup1_targetRxn.csv')
# annotation2$isTup1target = ifelse(colnames(rows.cor) %in% Tup1_rxn$Var1, 'Yes','No')
# colorList$isTup1target = c('#FD8D3CFF','#370335FF')
# names(colorList$isTup1target) = c('Yes','No')
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2[,c('isPredicted','isCorrelated','isTup1target','Tup1')],
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)

# plus Rap1
# Rap1_expression = read.csv('output/Rap1_expression.csv',row.names = 1)
# hist(cor(t(Rap1_expression),t(heatTbl)))
# corV = t(cor(t(Rap1_expression),t(heatTbl)))
# annotation2[rownames(corV),'Rap1'] = corV[,1]
# Rap1_rxn = read.csv('output/Rap1_targetRxn.csv')
# annotation2$isRap1target = ifelse(colnames(rows.cor) %in% Rap1_rxn$Var1, 'Yes','No')
# colorList$isRap1target = c('#FD8D3CFF','#370335FF')
# names(colorList$isRap1target) = c('Yes','No')
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2,
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)

#############covered by following analysis########
# # check Gln3
# Gln3_rxn = read.csv('output/Gln3_targetRxn.csv')
# annotation2$isGln3target = ifelse(colnames(rows.cor) %in% Gln3_rxn$Var1, 'Yes','No')
# colorList$isGln3target = c('#FD8D3CFF','#370335FF')
# names(colorList$isGln3target) = c('Yes','No')
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2[,c('isPredicted','isCorrelated','isGln3target')],
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)
# 
# 
# # check Gcn4
# Gcn4_rxn = read.csv('output/Gcn4_targetRxn.csv')
# annotation2$isGcn4target = ifelse(colnames(rows.cor) %in% Gcn4_rxn$Var1, 'Yes','No')
# colorList$isGcn4target = c('#FD8D3CFF','#370335FF')
# names(colorList$isGcn4target) = c('Yes','No')
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2[,c('isPredicted','isCorrelated','isGcn4target')],
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)
# 
# # check Lys14
# Lys14_rxn = read.csv('output/Lys14_targetRxn.csv')
# annotation2$isLys14target = ifelse(colnames(rows.cor) %in% Lys14_rxn$Var1, 'Yes','No')
# colorList$isLys14target = c('#FD8D3CFF','#370335FF')
# names(colorList$isLys14target) = c('Yes','No')
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2[,c('isPredicted','isCorrelated','isLys14target')],
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)
# 
# 
# # check Leu3
# Leu3_rxn = read.csv('output/Leu3_targetRxn.csv')
# annotation2$isLeu3target = ifelse(colnames(rows.cor) %in% Leu3_rxn$Var1, 'Yes','No')
# colorList$isLeu3target = c('#FD8D3CFF','#370335FF')
# names(colorList$isLeu3target) = c('Yes','No')
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2[,c('isPredicted','isCorrelated','isLeu3target')],
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)


########### check canonical TF in metabolism########
# label the coexpression 
RNArxnExp = read.csv('output/RNA_rxn_expression_FC_table.csv',row.names = 1)
RNAallExp = read.csv('output/RNA_all_expression_FC_table.csv')
# expression matrix for target rxns
rxnExpMat = t(RNArxnExp[colnames(rows.cor),])
colnames(rxnExpMat) = colnames(rows.cor)
rownames(RNAallExp) = RNAallExp$YORF

# load the GRN 
GRN = read.csv('output/regulationMatrix_metabolicOnly_matrix_rxn.csv',header = T)
GRN_rxnNames = read.csv('output/regulationMatrix_metabolicOnly_targetRxnNames.csv',header = T)
colnames(GRN) = GRN_rxnNames$Var1
geneIDtbl = read.csv('./../input/YeastJoshua/yeastGeneID.csv')
name2ID = data.frame(name = character(),ID = character())
for (i in 1:nrow(geneIDtbl)){
  IDs = strsplit(geneIDtbl$Gene.designations[i],';')
  name2ID = rbind(name2ID,data.frame(name = IDs[[1]],ID = geneIDtbl$OLN[i]))
}
allTF = read.csv('output/regulationMatrix_metabolicOnly_TFnames.csv',header = T)
allTF$TFid = allTF$TFnames
# intersect the selected TF and assign rownames
allTF$TFid[allTF$TFnames %in% name2ID$name] = name2ID$ID[match(allTF$TFnames[allTF$TFnames %in% name2ID$name], name2ID$name)]
allTF$TFid[allTF$TFnames %in% name2ID$ID] = name2ID$ID[match(allTF$TFnames[allTF$TFnames %in% name2ID$ID], name2ID$ID)]
allTF$TFid[str_replace(toupper(str_replace(allTF$TFnames,'p$','')),'-','_') %in% str_replace(toupper(name2ID$name),'-','_')] = name2ID$ID[match(
  str_replace(toupper(str_replace(allTF$TFnames[str_replace(toupper(str_replace(allTF$TFnames,'p$','')),'-','_') %in% str_replace(toupper(name2ID$name),'-','_')],'p$','')),'-','_'), 
  str_replace(toupper(name2ID$name),'-','_'))]
rownames(GRN) = allTF$TFid

# case-study#1 gcn4 -- the most famous master regulator of metabolism
# cgn4
# GCN4 YEL009C
Gcn4_expression = t(RNAallExp['YEL009C',3:ncol(RNAallExp)])
hist(cor(Gcn4_expression,rxnExpMat))

corV = t(cor(Gcn4_expression, rxnExpMat))
annotation2[rownames(corV),'Gcn4_coexpression_r'] = corV[,1]
colorList$Gcn4_coexpression_r = c('#00FF00','#ff0000')
names(colorList$Gcn4_coexpression_r) = c('-1','1')

pval = numeric()
for (i in 1:ncol(rxnExpMat)){
  if (!all(is.na(rxnExpMat[,i]))){
   tmp  = cor.test(Gcn4_expression, rxnExpMat[,i])
   pval[i] = tmp$p.value
  }else{
    pval[i] = NA
  }
}
fdrVal = p.adjust(pval,method = 'BH')
annotation2$isGcn4_coexpressed = ifelse(fdrVal < 0.05, 'Yes','No')
colorList$isGcn4_coexpressed = c('#FD8D3CFF','#370335FF')
names(colorList$isGcn4_coexpressed) = c('Yes','No')

Gcn4_rxn = read.csv('output/Gcn4_targetRxn.csv')
annotation2$isGcn4_target = ifelse(colnames(rows.cor) %in% Gcn4_rxn$Var1, 'Yes','No')
colorList$isGcn4_target = c('#FD8D3CFF','#370335FF')
names(colorList$isGcn4_target) = c('Yes','No')

pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,
         annotation_col = annotation2[,c('pathway_major','isCorrelated','isPredicted','isGcn4_target','isGcn4_coexpressed','Gcn4_coexpression_r')],
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off()
pdf('figures/GCN4_co-expression_coFlux.pdf',width = 24,height = 16)
pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,labels_col = labels_col,
         annotation_col = annotation2[,c('pathway_major','isCorrelated','isPredicted','isGcn4_target','isGcn4_coexpressed','Gcn4_coexpression_r')],
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 2, fontsize_col = 4,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off()
#gcn2 YDR283C ==> some pos corr in sugar metabolism, and strong neg corr in aa and puring pathways, alike gcn4!!!
#bas1 (YKR099W) ==> good positive coexpression with amino acid metabolism 
#bas2(pho2) YDL106C ==> very strong positive coexpression with some unlabeled rxns
#his4 YCL030C ==> very good and strong positive correlation with many amino acid pathways
#leu3 YLR451W ==> neg corr in puring and aa 
#lys14 YDR034C ==> lysine target but does not coexpress
#met4 YNL103W ==> many corr; pos with FA; both pos and neg with AA and puring 
#met28 YIR017C ==> few corr
#gln3  YER040W  ==> most amino acid pos cor; purine and lys not corr 

# selected examples: 

# indicating some relation between coexpression and correlatiojn
#gln3  YER040W  ==> very good and strong positive correlation with many amino acid pathways; purine and lys not corr 
#his4 YCL030C ==> very good and strong positive correlation with many amino acid pathways, especially histidine
#bas1 (YKR099W) ==> good positive coexpression with amino acid metabolism 
#leu3 YLR451W ==> neg corr in puring and aa, similar to gcn4 

# indicating mRNA level co-expression may not explain things
#lys14 YDR034C ==> lysine target but does not coexpress


targetTF = 'YDR034C'
outName = 'lys14'

targetTF_expression = t(RNAallExp[targetTF,3:ncol(RNAallExp)])
hist(cor(targetTF_expression,rxnExpMat))
corV = t(cor(targetTF_expression, rxnExpMat))
annotation2[rownames(corV),'TF_coexpression_r'] = corV[,1]
colorList$TF_coexpression_r = c('#00FF00','#ff0000')
names(colorList$TF_coexpression_r) = c('-1','1')
pval = numeric()
for (i in 1:ncol(rxnExpMat)){
  if (!all(is.na(rxnExpMat[,i]))){
    tmp  = cor.test(targetTF_expression, rxnExpMat[,i])
    pval[i] = tmp$p.value
  }else{
    pval[i] = NA
  }
}
fdrVal = p.adjust(pval,method = 'BH')
annotation2$is_coexpressed = ifelse(fdrVal < 0.05, 'Yes','No')
colorList$is_coexpressed = c('#FD8D3CFF','#370335FF')
names(colorList$is_coexpressed) = c('Yes','No')
TF_expression_rxn = colnames(GRN)[GRN[targetTF,] == 1]
annotation2$is_target = ifelse(colnames(rows.cor) %in% TF_expression_rxn, 'Yes','No')
colorList$is_target = c('#FD8D3CFF','#370335FF')
names(colorList$is_target) = c('Yes','No')
dev.off()
pdf(paste('figures/',outName,'_co-expression_coFlux.pdf',sep = ''),width = 24,height = 16)
pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,labels_col = labels_col,
         annotation_col = annotation2[,c('pathway_major','isPredicted','isCorrelated','is_target','is_coexpressed','TF_coexpression_r')],
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 2, fontsize_col = 4,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off()

##############didn find interesting target(too many)#######
# # search for the TF that is pos corr with puring - use greedy TF list
# allTF = read.csv('./../input/YeastJoshua/TFnames_fromRegulationMatrix.csv',header = F)
# geneIDtbl = read.csv('./../input/YeastJoshua/yeastGeneID.csv')
# name2ID = data.frame(name = character(),ID = character())
# for (i in 1:nrow(geneIDtbl)){
#   IDs = strsplit(geneIDtbl$Gene.designations[i],';')
#   name2ID = rbind(name2ID,data.frame(name = IDs[[1]],ID = geneIDtbl$OLN[i]))
# }
# allTF_ID1 = name2ID$ID[name2ID$name %in% allTF$V1]
# allTF_ID2 = name2ID$ID[name2ID$ID %in% allTF$V1]
# allTF_ID3 = name2ID$ID[str_replace(toupper(name2ID$name),'-','_') %in% str_replace(toupper(str_replace(allTF$V1,'p$','')),'-','_')]
# allTF_ID = unique(c(allTF_ID1,allTF_ID2,allTF_ID3))
# allTF_ID = intersect(allTF_ID, rownames(RNAallExp))
# puringrxns = rownames(annotation)[annotation$pathway_major == 'Purine metabolism']
# # check for correlation
# aveCorr = numeric()
# for (i in 1:length(allTF_ID)){
#   targetTF = allTF_ID[i]
#   targetTF_expression = t(RNAallExp[targetTF,3:ncol(RNAallExp)])
#   corV = t(cor(targetTF_expression,rxnExpMat[,puringrxns]))
#   aveCorr = c(aveCorr,mean(corV))
# }
# names(aveCorr) = allTF_ID
# aveCorr = sort(aveCorr,decreasing = T)
# aveCorr[1:10]
##############didn find interesting target(too many)#######
# search for the TF that is pos corr with puring - use TF that know to bind DNA with a motif
allTF = read.csv('./../input/YeastJoshua/TFnames_fromyetfasco.csv',header = T)
allTF_ID = intersect(allTF$SysName, rownames(RNAallExp))
puringrxns = rownames(annotation)[annotation$pathway_major == 'Purine metabolism']
# check for correlation
aveCorr = numeric()
for (i in 1:length(allTF_ID)){
  targetTF = allTF_ID[i]
  targetTF_expression = t(RNAallExp[targetTF,3:ncol(RNAallExp)])
  corV = t(cor(targetTF_expression,rxnExpMat[,puringrxns]))
  aveCorr = c(aveCorr,mean(corV))
}
names(aveCorr) = allTF_ID
aveCorr = sort(aveCorr,decreasing = T)
aveCorr[1:10]

# plot an interesting TF
targetTF = 'YDR451C'

targetTF_expression = t(RNAallExp[targetTF,3:ncol(RNAallExp)])
hist(cor(targetTF_expression,rxnExpMat))
corV = t(cor(targetTF_expression, rxnExpMat))
annotation2[rownames(corV),'TF_coexpression_r'] = corV[,1]
colorList$TF_coexpression_r = c('#00FF00','#ff0000')
names(colorList$TF_coexpression_r) = c('-1','1')
pval = numeric()
for (i in 1:ncol(rxnExpMat)){
  if (!all(is.na(rxnExpMat[,i]))){
    tmp  = cor.test(targetTF_expression, rxnExpMat[,i])
    pval[i] = tmp$p.value
  }else{
    pval[i] = NA
  }
}
fdrVal = p.adjust(pval,method = 'BH')
annotation2$is_coexpressed = ifelse(fdrVal < 0.05, 'Yes','No')
colorList$is_coexpressed = c('#FD8D3CFF','#370335FF')
names(colorList$is_coexpressed) = c('Yes','No')
# TF_expression_rxn = read.csv('output/Gcn4_targetRxn.csv')
# annotation2$isGcn4_target = ifelse(colnames(rows.cor) %in% Gcn4_rxn$Var1, 'Yes','No')
# colorList$isGcn4_target = c('#FD8D3CFF','#370335FF')
# names(colorList$isGcn4_target) = c('Yes','No')
pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,
         annotation_col = annotation2[,c('isPredicted','isCorrelated','is_coexpressed','TF_coexpression_r')],
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
##############search for the TF that is pos corr with purine and - use greedy TF list - known metabolic regulator#######
# search for the TF that is pos corr with puring - use greedy TF list - known metabolic regulator
allTF_ID = allTF$TFid
allTF_ID = intersect(allTF_ID, rownames(RNAallExp))
puringrxns = rownames(annotation)[annotation$pathway_major == 'Purine metabolism']

# intersect the selected TF and assign rownames
# check for correlation
aveCorr = numeric()
for (i in 1:length(allTF_ID)){
  targetTF = allTF_ID[i]
  targetTF_expression = t(RNAallExp[targetTF,3:ncol(RNAallExp)])
  corV = t(cor(targetTF_expression,rxnExpMat[,puringrxns]))
  # weight the cor r by the GRN 
  #corV = corV * exp(GRN[targetTF,rownames(corV)])
  aveCorr = c(aveCorr,sum(as.numeric(corV)))
}
names(aveCorr) = allTF_ID
aveCorr = sort(aveCorr,decreasing = T)
aveCorr[1:10]


# plot an interesting TF ==> there is no single magic TF that specifically targets and coexpress with purine. it remains an open question
# by weighted method
# YLR403W SFP1
# YIL131C FKH1: good correlation and also is target 
# YHR206W SKN7

# by unweighted
# YDR451C YHP1 very strong correlation, but not target genes
# YIL131C FKH1 
# YJR147W HMS2: paralog of SKN7, good correlation but not target 
# YGL162W correlated with most aa and purine genes

targetTF = 'YIL131C'
outName = 'FKH1'

targetTF_expression = t(RNAallExp[targetTF,3:ncol(RNAallExp)])
hist(cor(targetTF_expression,rxnExpMat))
corV = t(cor(targetTF_expression, rxnExpMat))
annotation2[rownames(corV),'TF_coexpression_r'] = corV[,1]
colorList$TF_coexpression_r = c('#00FF00','#ff0000')
names(colorList$TF_coexpression_r) = c('-1','1')
pval = numeric()
for (i in 1:ncol(rxnExpMat)){
  if (!all(is.na(rxnExpMat[,i]))){
    tmp  = cor.test(targetTF_expression, rxnExpMat[,i])
    pval[i] = tmp$p.value
  }else{
    pval[i] = NA
  }
}
fdrVal = p.adjust(pval,method = 'BH')
annotation2$is_coexpressed = ifelse(fdrVal < 0.05, 'Yes','No')
colorList$is_coexpressed = c('#FD8D3CFF','#370335FF')
names(colorList$is_coexpressed) = c('Yes','No')

TF_expression_rxn = colnames(GRN)[GRN[targetTF,] == 1]
annotation2$is_target = ifelse(colnames(rows.cor) %in% TF_expression_rxn, 'Yes','No')
colorList$is_target = c('#FD8D3CFF','#370335FF')
names(colorList$is_target) = c('Yes','No')
dev.off()
pdf(paste('figures/',outName,'_co-expression_coFlux.pdf',sep = ''),width = 24,height = 16)
pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("pathway_major"),drop = F],
         annotation_colors = colorList,labels_col = labels_col,
         annotation_col = annotation2[,c("pathway_major",'isPredicted','isCorrelated','is_target','is_coexpressed','TF_coexpression_r')],
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 2, fontsize_col = 4,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off()

########
# cross-correlation 
# it reveals that there are reactions whose expression does not correlate with its own flux but with some other flux module; this indicates
# that the reaction's gene is coexpressed with the genes of the other flux module; but the flux is not following the pattern 

r_mat = read.csv('output/crossCorrelation_rMat.csv',row.names = 1) # the pearson r of expression (row) -> flux (col)
colnames(r_mat) = rownames(r_mat)
#rows.cor <- cor(t(r_mat), use = "pairwise.complete.obs", method = "pearson")
#rows.dist <- dist(t(r_mat))
#library(biclust)
#rows.dist <- biclust(x = t(r_mat))

pheatmap(r_mat, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                    "RdYlBu")))(2000),
         #clustering_distance_rows = 'correlation',
         #clustering_distance_cols = 'correlation',
         annotation_row = annotation2[,c("pathway_major",'isPredicted','isCorrelated'),drop = F],
         annotation_col = annotation2[,c("pathway_major",'isPredicted','isCorrelated'),drop = F],
         annotation_colors = colorList,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE)

dev.off()
labels_col2 = paste(colnames(r_mat),annotation[colnames(r_mat),'pathway'])
labels_row2 = paste(rownames(r_mat),annotation[rownames(r_mat),'pathway'])
pdf('figures/cross-flux-expression-correlation.pdf',width = 28,height = 20)
pheatmap(r_mat, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                   "RdYlBu")))(2000),
         #clustering_distance_rows = 'correlation',
         #clustering_distance_cols = 'correlation',
         annotation_row = annotation2[,c("pathway_major",'isPredicted','isCorrelated'),drop = F],
         annotation_col = annotation2[,c("pathway_major",'isPredicted','isCorrelated'),drop = F],
         labels_col = labels_col2, labels_row = labels_row2,
         annotation_colors = colorList,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 7, fontsize_col = 7,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE,
         fontsize_number = 5)
dev.off()

# by p-value
p_mat = read.csv('output/crossCorrelation_pMat.csv',row.names = 1) # the pearson r of expression (row) -> flux (col)
colnames(p_mat) = rownames(p_mat)
p_mat = -log10(p_mat) * sign(r_mat)
pheatmap(p_mat, breaks = seq(-5,5,0.1),color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101),
         #clustering_distance_rows = rows.dist,
         #clustering_distance_cols = rows.dist,
         annotation_row = annotation2[,c("pathway_major",'isPredicted','isCorrelated'),drop = F],
         annotation_col = annotation2[,c("pathway_major",'isPredicted','isCorrelated'),drop = F],
         annotation_colors = colorList,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 5, fontsize_col = 8,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = TRUE,cluster_cols = TRUE)

# pairwise correlation of pattern for expression 
# Pairwise correlation between rows (genes)
# rows.cor <- cor(t(r_mat), use = "pairwise.complete.obs", method = "pearson")
# pheatmap(rows.cor, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                       "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation[,c("pathway_major"),drop = F],
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 5, fontsize_col = 8,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)

