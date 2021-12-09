# this script analyzed the pathway-level coexpression and pathway-level flux-expression correlation 
library(pheatmap)
library(RColorBrewer)
library(matrixStats)

# load data
boundary = data.frame(Var1 = c(c('base2-boundary6','expression only'), seq(0,40,0.5)))
FPAtbl = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
PCCs = read.csv('output/PCC_titration_all.csv',row.names = 1)
colnames(PCCs) = boundary$Var1
FPAtbl$maxPCC = rowMaxs(as.matrix(PCCs[rownames(FPAtbl),3:ncol(PCCs)]))
# load pathway annotation data
library(stringr)
library(xlsx)
annotation = read.xlsx('./pathway_annotations.xlsx','more_info',header = T)
rownames(annotation) = annotation$rxn

#correlation between coexpression levels and flux-expression concordance levels
ExpMat = read.csv('output/supp1D_normalizedExpression.csv',row.names = 1)
pathways = setdiff(unique(annotation$connected_pathways),'NA')
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
  subExp = ExpMat[rownames(ExpMat) %in% rxnset,]
  subExp = subExp[!duplicated(subExp),] # we only analyze the coexpression between valid (distinct) genes (GPR)
  
  if (nrow(subExp)>1){# only analyze the valid coexpression (at least two measured rxns [diff genes])
    tmp0 = cor(t(subExp), use = "pairwise.complete.obs", method = "pearson")
    tmp = c()# we avoid calculate the same PCC twice
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
excl = c('Leucine and valine biosynthesis')
cmpTbl_fit = cmpTbl[!(cmpTbl$pathways %in% excl), ]
cmpTbl_excl = cmpTbl[(cmpTbl$pathways %in% excl), ]

# PLOT
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


# next, we analyze the pathway level improvement from integration
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
model <- lm(y ~ poly(x,5)) # a simple fitting of trend line
seq1 = seq(-0.3,1, 0.001)
lines(seq1,predict(model,data.frame(x = seq1)),col='grey',lwd=1)
abline(h = 0,lty = 2)
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
box(lwd=1)
dev.off()
