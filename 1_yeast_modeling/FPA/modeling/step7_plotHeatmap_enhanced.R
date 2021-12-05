library(pheatmap)
library(RColorBrewer)
# the annotation of reactions
# the manual connected pathway was used; pathways with reactions less than 3 were grouped to others; isolated transporters are labeled out (inner and exchange)
mat = read.csv('output/relCorr_heatmapTbl_realDist.csv')
sig = read.csv('output/heatmapTbl_sigLabel_realDist.csv')
boundary = read.csv('output/heatmapTbl_boundaries_realDist.csv')

heatTbl = mat[,2:ncol(mat)]
colnames(heatTbl) = boundary$Var1[1:(ncol(mat)-1)]
rownames(heatTbl) = mat$Row
sigTbl = ifelse(sig[,2:ncol(mat)]==1, '+','')
dev.off()
pdf('figures/boundary_heatmap_realDist.pdf',width = 30,height = 17.5)
# we do row-wise left nearest impute: the prediction should be same to the left-most one if no new reaction can be seen as the distance boundary increase
for (i in 1:nrow(heatTbl)){
  for (j in 1:(ncol(heatTbl))){
    if (is.na(heatTbl[i,j])){
      p1 = heatTbl[i,1:j-1]
      p1 = p1[!is.na(p1)]
      s1 = sigTbl[i,1:j-1]
      s1 = s1[!is.na(p1)]
      heatTbl[i,j] = p1[length(p1)]#mean(c(p1[length(p1)],p2[1]))
      sigTbl[i,j] = s1[length(s1)]#mean(c(p1[length(p1)],p2[1]))
    }
  }
}

# reorder according to annotation table 
library(xlsx)
annTbl = read.xlsx('./pathway_annotations.xlsx','all',header = T)
rownames(annTbl) = annTbl$rxn
annTbl = annTbl[rownames(annTbl) %in% rownames(heatTbl),]
       
annotationList = c('Purine metabolism','Phenylalanine, tyrosine and tryptophan biosynthesis','Histidine metabolism','Lysine biosynthesis',
                   'Arginine biosynthesis','Threonine, methionine and cysteine synthesis','Proline biosynthesis','Pyrimidine metabolism',
                   'Glycolysis','TCA cycle', 'Mannan synthesis',
                   'UDP-D-glucose metabolism','Fatty acid biosynthesis',
                   'Transporter [inner]','Transporter [exchange]','Others')
colors = read.csv('figures/simpsons_color.csv',header = F)
colorList2 = colors$V1[2:(1+length(annotationList))]
names(colorList2) = annotationList
colorList2 = list(functional_annotation = colorList2)
colorList2$functional_annotation = colorList2$functional_annotation[unique(annTbl[,c("functional_annotation")])]
# assign annotations
rownames(sigTbl) = rownames(heatTbl)
heatTbl = heatTbl[annTbl$rxn,]
sigTbl = sigTbl[annTbl$rxn,]
ann2 = annTbl[rownames(annTbl),"connected_pathway"]
ann2[ann2 == 'NA'] = annTbl$pathway[ann2 == 'NA']
labels_row0 = paste(rownames(annTbl),ann2)

pheatmap(heatTbl[,1:23], breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                      "RdYlBu")))(1000),
         annotation_row = annTbl[,c("functional_annotation"),drop = F],
         annotation_colors = colorList2, labels_row = labels_row0,
         gaps_col = c(1,2),
         border_color = 'black',
         fontsize = 8,fontsize_row = 6, fontsize_col = 8,
         cellwidth = 12, cellheight = 6,
         cluster_rows = FALSE,cluster_cols = FALSE
         ,display_numbers = sigTbl[,1:23], number_color = 'Black',fontsize_number = 6
)

dev.off()

# plot the subset showcase
library(parallel)
smooth.expr <- function(mat, x, x.pred) {
  n <- length(x.pred)
  nd <- data.frame(x=x.pred)
  tmp <- mclapply(1:nrow(mat), function(i) {
    fit <- loess(mat[i, ] ~ x, span=0.1, degree=2)
    predict(fit, newdata=nd, se=FALSE)
  })
  expr.fit <- t(matrix(unlist(tmp), n))
  rownames(expr.fit) <- rownames(mat)
  colnames(expr.fit) <- nd$x
  return(expr.fit)
}
library(lsa)
# smooth to make the trend of change more visually clear
annTbl = read.xlsx('./pathway_annotations.xlsx','selected',header = T)
rownames(annTbl) = annTbl$rxn
annTbl = annTbl[rownames(annTbl) %in% rownames(heatTbl),]
heatTbl = heatTbl[annTbl$rxn,]
sigTbl = sigTbl[annTbl$rxn,]
labels_row0 = paste(rownames(annTbl),annTbl[rownames(annTbl),"functional_annotation"])
tg_fit <- smooth.expr(as.matrix(heatTbl[,3:23]),as.numeric(colnames(heatTbl)[3:23]),seq(0,20,0.01))
gaps = c(0)
for (i in 1:(nrow(annTbl)-1)){
  if(annTbl$functional_annotation[i] != annTbl$functional_annotation[i+1]){
    gaps = c(gaps, i)
  }
}
colorList3 = colorList2
colorList3$functional_annotation = colorList3$functional_annotation[unique(annTbl[,c("functional_annotation")])]
dev.off()
png('figures/boundary_heatmap_realDist_smooth.png',width = 1200,height = 800)
pheatmap(tg_fit, breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                           "RdYlBu")))(1000),
         annotation_row = annTbl[,c("functional_annotation"),drop = F],
         annotation_colors = colorList3, #labels_row = labels_row0,
         show_colnames = F, show_rownames = F,
         border_color = 0,
         gaps_row = gaps,
         cluster_rows = FALSE,cluster_cols = FALSE
)
dev.off()

