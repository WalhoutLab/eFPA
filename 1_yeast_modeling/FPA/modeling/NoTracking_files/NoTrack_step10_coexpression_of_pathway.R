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
cmpTbl$max_coexpression = NA
cmpTbl$ave_FPA_corr = NA
cmpTbl$median_FPA_corr = NA
cmpTbl$ave_corr = NA
cmpTbl$median_corr = NA
cmpTbl$max_corr = NA
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
    cmpTbl$max_coexpression[i] =  max(tmp,na.rm = T)
    cmpTbl$ave_FPA_corr[i] = mean(FPAtbl[rxnset,'PCC_by_default_FPA'],na.rm = T)
    cmpTbl$median_FPA_corr[i] = median(FPAtbl[rxnset,'PCC_by_default_FPA'],na.rm = T)
    cmpTbl$ave_corr[i] = mean(FPAtbl[rxnset,"PCC_by_target_expression"],na.rm = T)
    cmpTbl$median_corr[i] = median(FPAtbl[rxnset,'PCC_by_target_expression'],na.rm = T)
    cmpTbl$max_corr[i] = max(FPAtbl[rxnset,'PCC_by_target_expression'],na.rm = T)
    cmpTbl$ave_max_corr[i] = mean(FPAtbl[rxnset,"maxPCC"],na.rm = T)
    cmpTbl$median_max_corr[i] = median(FPAtbl[rxnset,'maxPCC'],na.rm = T)
  }
}
cmpTbl = cmpTbl[!is.na(cmpTbl$median_coexpression),]
excl = c('Leucine and valine biosynthesis')
cmpTbl_fit = cmpTbl[!(cmpTbl$pathways %in% excl), ]
cmpTbl_excl = cmpTbl[(cmpTbl$pathways %in% excl), ]

# PLOT
pdf('figures/NoTrack_median-concordance-max_concordance.pdf',width = 4,height = 4.3)
fit = lm(max_corr ~median_corr, cmpTbl_fit)
plot(cmpTbl_fit$median_corr, cmpTbl_fit$max_corr, ylim = c(-0.3,1), xlim = c(-0.7,1),
     xlab = 'Pathway-level medain flux-expression correlation',
     ylab = 'Pathway-level max flux-expression correlation',
     lwd = 1)
text(cmpTbl_fit$median_corr, cmpTbl_fit$max_corr, labels=rank(cmpTbl_fit$median_coexpression), cex = 0.3)

abline(fit,lwd = 1,lty = 2)
points(cmpTbl_excl$median_corr, cmpTbl_excl$max_corr,col = 'red',pch = 1, lwd = 1)
a = summary(fit)
text(0.8, -0.1,paste('r = ',round(sqrt(a$r.squared),2),sep = ''))
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
dev.off()


# PLOT
pdf('figures/NoTrack_optimal_FPA-max_concordance.pdf',width = 4,height = 4.3)
fit = lm(max_corr ~median_max_corr, cmpTbl_fit)
plot(cmpTbl_fit$median_max_corr, cmpTbl_fit$max_corr, ylim = c(-0.3,1), xlim = c(-0.7,1),
     xlab = 'Pathway-level optimal FPA prediction correlation',
     ylab = 'Pathway-level max flux-expression correlation',
     lwd = 1)
text(cmpTbl_fit$median_max_corr, cmpTbl_fit$max_corr, labels=rank(cmpTbl_fit$median_coexpression), cex = 0.3)

abline(fit,lwd = 1,lty = 2)
points(cmpTbl_excl$median_max_corr, cmpTbl_excl$max_corr,col = 'red',pch = 1, lwd = 1)
a = summary(fit)
text(0.8, -0.1,paste('r = ',round(sqrt(a$r.squared),2),sep = ''))
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
dev.off()

# ==> Pyrimidine metabolism and Leucine and valine biosynthesis  are stong control point 

# (1) we analyze the pathway level improvement from integration (for default FPA)
# NOTE: the result is qualitatively similar to this if we calculate the delta of median PCC only for the 156 measured reactions
cmpTbl2 = cmpTbl
cmpTbl2 = cmpTbl2[order(cmpTbl2$median_coexpression),]
dev.off()
pdf('figures/NoTrack_co-expression-flux-accordance_deltaPCC.pdf',width = 4,height = 4.3)
y = cmpTbl2$median_FPA_corr - cmpTbl2$median_corr
x = cmpTbl2$median_coexpression
plot(x, y, ylim = c(-0.25,0.4), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level integration benefit (delta PCC)',
     lwd = 1,cex = 0.5,
     type = 'p', lty =2)
model <- lm(y ~ poly(x,3)) # a simple fitting of trend line
seq1 = seq(-0.3,1, 0.001)
lines(seq1,predict(model,data.frame(x = seq1)),col='grey',lwd=1)
abline(h = 0,lty = 2)
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
box(lwd=1)
dev.off()

# (2) we check the delta median vs. median delta
pathways = cmpTbl$pathways
cmpTbl$median_deltaPCC_max = NA
cmpTbl$median_deltaPCC_FPA = NA
for (i in 1:length(pathways)){
  rxnset = annotation$rxn[which(annotation$connected_pathways == pathways[i])]
  cmpTbl$median_deltaPCC_max[i] = median(FPAtbl[rxnset,'maxPCC'] - FPAtbl[rxnset,'PCC_by_target_expression'],na.rm = T)
  cmpTbl$median_deltaPCC_FPA[i] = median(FPAtbl[rxnset,'PCC_by_default_FPA'] - FPAtbl[rxnset,'PCC_by_target_expression'],na.rm = T)
}
cmpTbl2 = cmpTbl
cmpTbl2 = cmpTbl2[order(cmpTbl2$median_coexpression),]
dev.off()
pdf('figures/NoTrack_co-expression-flux-accordance_deltaPCC_median_delta_FPA.pdf',width = 4,height = 4.3)
y = cmpTbl2$median_deltaPCC_FPA
x = cmpTbl2$median_coexpression
plot(x, y, ylim = c(-0.2,0.5), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level integration benefit (delta PCC)',
     lwd = 1,cex = 0.5,
     type = 'p', lty =2)
model <- lm(y ~ poly(x,3)) # a simple fitting of trend line
seq1 = seq(-0.3,1, 0.001)
lines(seq1,predict(model,data.frame(x = seq1)),col='grey',lwd=1)
abline(h = 0,lty = 2)
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
box(lwd=1)
dev.off()
dev.off()
pdf('figures/NoTrack_co-expression-flux-accordance_deltaPCC_median_delta_max.pdf',width = 4,height = 4.3)
y = cmpTbl2$median_deltaPCC_max
x = cmpTbl2$median_coexpression
plot(x, y, ylim = c(-0.1,1.1), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level integration benefit (delta PCC)',
     lwd = 1,cex = 0.5,
     type = 'p', lty =2)
model <- lm(y ~ poly(x,2)) # a simple fitting of trend line
seq1 = seq(-0.3,1, 0.001)
lines(seq1,predict(model,data.frame(x = seq1)),col='grey',lwd=1)
abline(h = 0,lty = 2)
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
box(lwd=1)
dev.off()
# ==> we found that median delta PCC disminished the trend and benefit for default FPA, and also made the max FPA worse (but still interpretable and trend is clear)

# (3) directly show if coexpression in curated pathway is predictive of flux 
FluxTbl = read.csv('output/supp1B_normalizedFlux.csv',row.names = 1)
FPAtbl$coexpressionPCC = NA
for (i in 1:length(pathways)){
  pathways[i]
  rxnset = annotation$rxn[which(annotation$connected_pathways == pathways[i])]
  subExp = ExpMat[rownames(ExpMat) %in% rxnset,]
  subExp = subExp[!duplicated(subExp),] # we only analyze the coexpression between valid (distinct) genes (GPR)
  subExp = t(scale(t(subExp)))
  
  if (nrow(subExp)>1){# only analyze the valid coexpression (at least two measured rxns [diff genes])
    # calculate coexpression 
    #subExp.pca <- prcomp(subExp)
    #trend = subExp.pca$rotation[,1]
    trend = colMedians(subExp)
    # the trend can have sign flip (loading's sign is arbitury), we fix it
    # if (median(cor(as.numeric(trend), t(subExp))) < 0){
    #   trend = -trend
    # }
    for (rxn in rxnset){
      FPAtbl[rxn, 'coexpressionPCC'] = cor(abs(as.numeric(FluxTbl[rxn,])), as.numeric(trend))
    }
    cor( t(subExp))
    cor(as.numeric(trend), t(subExp))
    FPAtbl[rxnset, 'coexpressionPCC'] 
    FPAtbl[rxnset, "PCC_by_target_expression"]
    FPAtbl[rxnset, "PCC_by_default_FPA"]
    FPAtbl[rxnset, "maxPCC"]
  }
}

#correlation between coexpression levels and PC1-expression concordance levels
cmpTbl$median_PC1_corr = NA
for (i in 1:length(pathways)){
  rxnset = annotation$rxn[which(annotation$connected_pathways == pathways[i])]
  cmpTbl$median_PC1_corr[i] = median(FPAtbl[rxnset,'coexpressionPCC'],na.rm = T)
}
excl = c('Leucine and valine biosynthesis')
cmpTbl_fit = cmpTbl[!(cmpTbl$pathways %in% excl), ]
cmpTbl_excl = cmpTbl[(cmpTbl$pathways %in% excl), ]

# PLOT
#dev.off()
#pdf('figures/co-expression-flux-accordance.pdf',width = 4,height = 4.3)
fit = lm(median_PC1_corr ~median_coexpression, cmpTbl_fit)
plot(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_PC1_corr, ylim = c(-1,1), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level flux-expression correlation',
     lwd = 1)
text(cmpTbl_fit$median_coexpression, cmpTbl_fit$median_PC1_corr, labels=rank(cmpTbl_fit$median_coexpression), cex = 0.3)

abline(fit,lwd = 1,lty = 2)
points(cmpTbl_excl$median_coexpression, cmpTbl_excl$median_PC1_corr,col = 'red',pch = 1, lwd = 1)
a = summary(fit)
text(0.8, -0.1,paste('r = ',round(sqrt(a$r.squared),2),sep = ''))
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)

# investigate the integration benefit
y = cmpTbl$median_PC1_corr - cmpTbl$median_corr
x = cmpTbl$median_coexpression
plot(x, y, ylim = c(-0.3,0.4), xlim = c(-0.3,1),
     xlab = 'Pathway-level co-expression',
     ylab = 'Pathway-level integration benefit (delta PCC)',
     lwd = 1,cex = 0.5,
     type = 'p', lty =2)
model <- lm(y ~ poly(x,4)) # a simple fitting of trend line
seq1 = seq(-0.3,1, 0.001)
lines(seq1,predict(model,data.frame(x = seq1)),col='grey',lwd=1)
abline(h = 0,lty = 2)
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
box(lwd=1)

#dev.off()
#==> coexpression trend (by the medians) is predictive of flux and correlated with the level of
# coexpression. however, the benefit of (optimal bound) FPA does not always come from the coexpression
# infomation of the target pathway itself, could also be distal (the local FPA usually dependents on coexpression 
# of the target pathway). so, the idea of coexpression dicatating flux is possible right, but extracting the trend
# of coexpression may not always improve the prediction.






