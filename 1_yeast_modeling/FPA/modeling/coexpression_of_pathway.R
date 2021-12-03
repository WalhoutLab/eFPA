library(pheatmap)
library(RColorBrewer)
library(matrixStats)
mat = read.csv('output/relCorr_heatmapTbl_wtdDist.csv')
sig = read.csv('output/heatmapTbl_sigLabel.csv')
boundary = read.csv('output/heatmapTbl_boundaries.csv')

FPAtbl = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
sigTbl = ifelse(sig[,2:ncol(mat)]==1, '+','')

# relPCCtbl = mat[,2:ncol(mat)]
# colnames(relPCCtbl) = boundary$Var1[1:(ncol(mat)-1)]
# rownames(relPCCtbl) = mat$Row
# relPCCtbl$absPCC_default = FPAtbl[rownames(relPCCtbl),"PCC_by_default_FPA"]
# relPCCtbl$maxAbsPCC = relPCCtbl$absPCC_default / relPCCtbl$`base 2 - boundary 6`
PCCs = read.csv('output/PCC_titration_all.csv',row.names = 1)
colnames(PCCs) = boundary$Var1[1:(ncol(mat)-1)]
FPAtbl$maxPCC = rowMaxs(as.matrix(PCCs[rownames(FPAtbl),3:ncol(PCCs)]))

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

########the coexpression heatmap of proteins(rxns)########
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
# the coexpression of fluxes
FluxMat = read.csv('output/relativeFluxTable.csv',row.names = 1)
sigInModel = mat$Row
# Pairwise correlation between rows (genes)
rows.cor <- cor(t(abs(FluxMat)), use = "pairwise.complete.obs", method = "pearson")
# prepare the new annotation table 
colorList$isPredicted = c('#FD8D3CFF','#370335FF')
names(colorList$isPredicted) = c('Yes','No')
colorList$isCorrelated = c('#FD8D3CFF','#370335FF')
names(colorList$isCorrelated) = c('Yes','No')
annotation2 = data.frame(row.names = colnames(rows.cor), isPredicted = ifelse(colnames(rows.cor) %in% sigInModel,'Yes','No'),
                         isCorrelated = ifelse(colnames(rows.cor) %in% mat$Row[sigTbl[,boundary == 'expression only'] == '+'], 'Yes','No'))
annotation2$pathway_major = annotation[rownames(annotation2), 'pathway_major']
annotation2$connected_pathways = annotation[rownames(annotation2), 'connected_pathways']
annotation2 = annotation2[,c('pathway_major','isCorrelated','isPredicted','connected_pathways')]
annotation2$pearson_r = pearson_r[rownames(annotation2),1]
colorList$pearson_r = c('#00FF00','#ff0000')
names(colorList$pearson_r) = c('-1','1')
rows.cor.exp <- cor(t(ExpMat), use = "pairwise.complete.obs", method = "pearson")
labels_col2 = paste(colnames(rows.cor.exp),annotation[colnames(rows.cor.exp),'pathway'])
inlabels = annotation2[rownames(rows.cor.exp),'connected_pathways']
rows.cor.exp = rows.cor.exp[order(inlabels),order(inlabels)]
pheatmap(rows.cor.exp, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                      "RdYlBu")))(2000),
         # clustering_distance_rows = as.dist(1 - rows.cor),
         annotation_row = annotation[,c("connected_pathways"),drop = F],
         annotation_colors = colorList,
         annotation_col = annotation2[,c("connected_pathways",'pearson_r','isCorrelated','isPredicted')],labels_col = labels_col2,
         gaps_col = 1,
         # legend_breaks = seq2,legend_labels = legend_labels,
         border_color = 'black',
         fontsize = 8,fontsize_row = 3, fontsize_col = 7,
         #cellwidth = 8, cellheight = 4,
         cluster_rows = F,cluster_cols = F,
         fontsize_number = 5)

######### consider the pathway-level co-expression matrix ######
ExpMat = read.csv('output/relativeExpressionTable.csv',row.names = 1)
rows.cor.exp <- cor(t(ExpMat), use = "pairwise.complete.obs", method = "pearson")
pathways = setdiff(unique(annotation$connected_pathways),NA)
cmpTbl = data.frame(pathways)
cmpTbl$ave_coexpression = NA
cmpTbl$median_coexpression = NA
cmpTbl$ave_FPA_corr = NA
cmpTbl$median_FPA_corr = NA
cmpTbl$ave_corr = NA
cmpTbl$median_corr = NA
cmpTbl$ave_max_corr = NA
cmpTbl$median_max_corr = NA
for (i in 1:length(pathways)){
  rxnset = annotation$rxn[which(annotation$connected_pathways == pathways[i])]
  if (sum(rownames(rows.cor.exp) %in% rxnset)>1){
    tmp0 = rows.cor.exp[rownames(rows.cor.exp) %in% rxnset, 
                       colnames(rows.cor.exp) %in% rxnset]
    tmp = c()
    for (zz in 1:(ncol(tmp0)-1)){
      for (kk in (zz+1):nrow(tmp0)){
        tmp = c(tmp, tmp0[zz,kk])
      }
    }
    
    cmpTbl$ave_coexpression[i] =  mean(tmp,na.rm = T)
    cmpTbl$median_coexpression[i] =  median(tmp,na.rm = T)
    cmpTbl$ave_FPA_corr[i] = mean(FPAtbl[rxnset,'PCC_by_default_FPA'],na.rm = T)
    cmpTbl$median_FPA_corr[i] = median(FPAtbl[rxnset,'PCC_by_default_FPA'],na.rm = T)
    cmpTbl$ave_corr[i] = mean(FPAtbl[rxnset,"PCC_by_target_expression"],na.rm = T)
    cmpTbl$median_corr[i] = median(FPAtbl[rxnset,'PCC_by_target_expression'],na.rm = T)
    cmpTbl$ave_max_corr[i] = mean(FPAtbl[rxnset,"maxPCC"],na.rm = T)
    cmpTbl$median_max_corr[i] = median(FPAtbl[rxnset,'maxPCC'],na.rm = T)
  }
}
cmpTbl = cmpTbl[!is.na(cmpTbl$median_coexpression),]
excl = c('leucine and valine biosynthesis')
cmpTbl_fit = cmpTbl[!(cmpTbl$pathways %in% excl), ]
cmpTbl_excl = cmpTbl[(cmpTbl$pathways %in% excl), ]

# plot for expression only
fit = lm(median_corr ~median_coexpression, cmpTbl_fit)
plot(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_corr, ylim = c(-1,1), xlim = c(0,1),
     xlab = 'median of all pairwise PCC for a pathway',
     ylab = 'median of FPA PCCs for a pathway')
abline(fit)
points(cmpTbl_excl$median_coexpression, cmpTbl_excl$median_corr,col = 'red')
a = summary(fit)
text(0.2, -0.6,paste('r squared = ',round(a$r.squared,2),sep = ''))

# plot for FPA
fit = lm(median_FPA_corr ~median_coexpression, cmpTbl_fit)
plot(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_FPA_corr, ylim = c(-1,1), xlim = c(0,1),
     xlab = 'median of all pairwise PCC for a pathway',
     ylab = 'median of FPA PCCs for a pathway')
abline(fit)
points(cmpTbl_excl$median_coexpression, cmpTbl_excl$median_FPA_corr,col = 'red')
a = summary(fit)
text(0.2, -0.6,paste('r squared = ',round(a$r.squared,2),sep = ''))

# plot for best FPA
fit = lm(median_max_corr ~median_coexpression, cmpTbl_fit)
plot(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_max_corr, ylim = c(-1,1), xlim = c(0,1),
     xlab = 'median of all pairwise PCC for a pathway',
     ylab = 'median of FPA max PCCs for a pathway')
abline(fit)
points(cmpTbl_excl$median_coexpression, cmpTbl_excl$median_max_corr,col = 'red')
a = summary(fit)
text(0.2, -0.6,paste('r squared = ',round(a$r.squared,2),sep = ''))


# all together
dev.off()
pdf('figures/co-expression-flux-accordance.pdf',width = 7,height = 7)
fit = lm(median_corr ~median_coexpression, cmpTbl_fit)
plot(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_corr, ylim = c(-1,1), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level flux-expression correlation',
     lwd = 2)
arrows(cmpTbl$median_coexpression, cmpTbl$median_FPA_corr,
       cmpTbl$median_coexpression, cmpTbl$median_max_corr, col = '#77AC30',
       length = 0.1)
arrows(cmpTbl$median_coexpression, cmpTbl$median_corr,
       cmpTbl$median_coexpression, cmpTbl$median_FPA_corr, col = '#0072BD',
       length = 0.1)
abline(fit,lwd = 2,lty = 2)
points(cmpTbl_excl$median_coexpression, cmpTbl_excl$median_corr,col = 'red',pch = 19, lwd = 2)
a = summary(fit)
# plot for FPA
points(cmpTbl$median_coexpression, cmpTbl$median_FPA_corr, col = '#0072BD', lwd = 2)
points(cmpTbl$median_coexpression, cmpTbl$median_max_corr, col = '#77AC30', lwd = 2)
legend(x = -0.35,y = 1.05,legend = c('expression','local integration','optimal boundary integration'),
       pch = 1,bty = 'n',lty=0,col = c('black','#0072BD','#77AC30'), cex=1,lwd = 2)
text(0.8, -0.1,paste('r = ',round(sqrt(a$r.squared),2),sep = ''))
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
box(lwd=2)
dev.off()


# expression only plot (publication)

# all together
dev.off()
pdf('figures/co-expression-flux-accordance.pdf',width = 4,height = 4.3)
fit = lm(median_corr ~median_coexpression, cmpTbl_fit)
plot(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_corr, ylim = c(-1,1), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level flux-expression correlation',
     lwd = 1)
text(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_corr, labels=rank(cmpTbl_fit$median_coexpression), cex = 0.3)

abline(fit,lwd = 1,lty = 2)
points(cmpTbl_excl$median_coexpression, cmpTbl_excl$median_corr,col = 'red',pch = 1, lwd = 1)
a = summary(fit)
text(0.8, -0.1,paste('r = ',round(sqrt(a$r.squared),2),sep = ''))
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
dev.off()
cmpTbl_fit$numberInFig = rank(cmpTbl_fit$median_coexpression)
pathwayCounts = data.frame(table(annotation$connected_pathways))
pathwayCounts$valid = 0
for (i in 1:nrow(pathwayCounts)){
  rxnset = annotation$rxn[which(annotation$connected_pathways == pathwayCounts$Var1[i])]
  pathwayCounts$valid[i]=sum(rownames(rows.cor.exp) %in% rxnset)
}


# pathway level improvement
cmpTbl2 = cmpTbl
cmpTbl2 = cmpTbl2[order(cmpTbl2$median_coexpression),]
dev.off()
pdf('figures/co-expression-flux-accordance_deltaPCC.pdf',width = 4,height = 4.3)
y = cmpTbl2$median_max_corr - cmpTbl2$median_corr
x = cmpTbl2$median_coexpression
plot(x, y, ylim = c(-0.1,1.05), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level integration benefit (delta PCC)',
     lwd = 1,cex = 0.5,
     type = 'p', lty =2)
model <- lm(y ~ poly(x,5))
seq1 = seq(-0.3,1, 0.001)
lines(seq1,predict(model,data.frame(x = seq1)),col='grey',lwd=1)
abline(h = 0,lty = 2)
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
box(lwd=1)
dev.off()




# dev.off()
# pdf('figures/co-expression.pdf',width = 24,height = 20)
# pheatmap(rows.cor.exp, breaks = seq(-1,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                                                           "RdYlBu")))(2000),
#          # clustering_distance_rows = as.dist(1 - rows.cor),
#          annotation_row = annotation[,c("pathway_major"),drop = F],
#          annotation_colors = colorList,
#          annotation_col = annotation2[,c("pathway_major",'pearson_r','isCorrelated','isPredicted')],labels_col = labels_col2,
#          gaps_col = 1,
#          # legend_breaks = seq2,legend_labels = legend_labels,
#          border_color = 'black',
#          fontsize = 8,fontsize_row = 3, fontsize_col = 7,
#          #cellwidth = 8, cellheight = 4,
#          cluster_rows = TRUE,cluster_cols = TRUE,
#          fontsize_number = 5)
# dev.off()
# 


