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

#correlation between coflux levels and flux-expression concordance levels
ExpMat = read.csv('output/supp1D_normalizedExpression.csv',row.names = 1)
fluxMat = read.csv('output/supp1B_normalizedFlux.csv',row.names = 1)
pathways = setdiff(unique(annotation$connected_pathways),'NA')
cmpTbl = data.frame(pathways)
cmpTbl$ave_coexpression = NA
cmpTbl$median_coexpression = NA
cmpTbl$ave_coflux = NA
cmpTbl$median_coflux = NA
cmpTbl$ave_FPA_corr = NA
cmpTbl$median_FPA_corr = NA
cmpTbl$ave_corr = NA
cmpTbl$median_corr = NA
cmpTbl$ave_max_corr = NA
cmpTbl$median_max_corr = NA
usedRxns = c()
for (i in 1:length(pathways)){
  rxnset = annotation$rxn[which(annotation$connected_pathways == pathways[i])]
  subExp = ExpMat[rownames(ExpMat) %in% rxnset,]
  subExp = subExp[!duplicated(subExp),] # we only analyze the coexpression between valid (distinct) genes (GPR)
  
  subFlux = abs(fluxMat[rownames(fluxMat) %in% rxnset,])

  if (nrow(subExp)>1){# only analyze the valid coexpression (at least two measured rxns [diff genes])
    tmp0 = cor(t(subExp), use = "pairwise.complete.obs", method = "pearson")
    tmp = c()# we avoid calculate the same PCC twice
    for (zz in 1:(ncol(tmp0)-1)){
      for (kk in (zz+1):nrow(tmp0)){
        tmp = c(tmp, tmp0[zz,kk])
      }
    }
    usedRxns = union(usedRxns, rownames(subExp))
    cmpTbl$ave_coexpression[i] =  mean(tmp,na.rm = T)
    cmpTbl$median_coexpression[i] =  median(tmp,na.rm = T)
    cmpTbl$ave_FPA_corr[i] = mean(FPAtbl[rxnset,'PCC_by_default_FPA'],na.rm = T)
    cmpTbl$median_FPA_corr[i] = median(FPAtbl[rxnset,'PCC_by_default_FPA'],na.rm = T)
    cmpTbl$ave_corr[i] = mean(FPAtbl[rxnset,"PCC_by_target_expression"],na.rm = T)
    cmpTbl$median_corr[i] = median(FPAtbl[rxnset,'PCC_by_target_expression'],na.rm = T)
    cmpTbl$ave_max_corr[i] = mean(FPAtbl[rxnset,"maxPCC"],na.rm = T)
    cmpTbl$median_max_corr[i] = median(FPAtbl[rxnset,'maxPCC'],na.rm = T)
  }
  
  if (nrow(subFlux)>1){# only analyze the valid coexpression (at least two measured rxns [diff genes])
    tmp0 = cor(t(subFlux), use = "pairwise.complete.obs", method = "pearson")
    tmp = c()# we avoid calculate the same PCC twice
    for (zz in 1:(ncol(tmp0)-1)){
      for (kk in (zz+1):nrow(tmp0)){
        tmp = c(tmp, tmp0[zz,kk])
      }
    }
    cmpTbl$ave_coflux[i] =  mean(tmp,na.rm = T)
    cmpTbl$median_coflux[i] =  median(tmp,na.rm = T)
  }
}
#write.csv(x = usedRxns,'output/rxns_used_for_coexp_analysis.csv')
cmpTbl = cmpTbl[!is.na(cmpTbl$median_coexpression),]
excl = c('Leucine and valine biosynthesis')
cmpTbl_fit = cmpTbl[!(cmpTbl$pathways %in% excl), ]
cmpTbl_excl = cmpTbl[(cmpTbl$pathways %in% excl), ]

# PLOT
dev.off()
pdf('figures/co-flux-flux-accordance.pdf',width = 3.3,height = 3.655)
fit = lm(median_corr ~median_coflux, cmpTbl_fit)
plot(cmpTbl_fit$median_coflux, cmpTbl_fit$median_corr, ylim = c(-0.7,1), xlim = c(0,1),
     xlab = 'Pathway-level co-flux',
     ylab = 'Pathway-level flux-expression correlation',
     lwd = 1)
text(cmpTbl_fit$median_coflux, cmpTbl_fit$median_corr, labels=rank(cmpTbl_fit$median_coexpression), cex = 0.3)
abline(fit,lwd = 1,lty = 2)
points(cmpTbl_excl$median_coflux, cmpTbl_excl$median_corr,col = 'red',pch = 1, lwd = 1)
a = summary(fit)
text(0.8, -0.1,paste('r = ',round(sqrt(a$r.squared),2),sep = ''))
axis(side = 1, lwd = 1)
axis(side = 2, lwd = 1)
dev.off()

