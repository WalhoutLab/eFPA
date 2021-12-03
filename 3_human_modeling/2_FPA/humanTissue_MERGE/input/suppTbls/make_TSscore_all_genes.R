# the original paper (Cell 2021) only provided the RNA TS score for genes quantified by both proteomics and 
# RNA-seq. We calculated the TS score for all genes detected in RNA-seq for downstream analysis
# we found a slight differene numerically when we tried to reproduce the TS score calculation in the original 
# paper. To keep minimal inconsistency of the data, we only use in-house calculated TS score for the genes that 
# the original paper didnt provide a TS score on
# In addition, to simplify our analysis, we didnt filter the new TS scores as the original paper did. 

source('AdaReg.R')

data = read.csv('rnaAllSamples_log2TPM.csv')

zMat = matrix(NA,nrow = nrow(data),ncol = ncol(data)-2)
rownames(zMat) = data$gene.id.full
colnames(zMat) = colnames(data)[3:ncol(data)]
fittingRes = data.frame(geneID = data$gene.id.full)
fittingRes$popMean = NA
fittingRes$popSD = NA
for (i in 1:nrow(data)){
  y = as.numeric(data[i,3:ncol(data)])
  out = AdaReg(model.matrix(~1,data = as.data.frame(y)), y)
  zr1 = (y-out$beta.rob.fit)/sqrt(out$var.sig.gp.fit)
  zMat[i,] = zr1
  fittingRes$popMean[i] = out$beta.rob.fit
  fittingRes$popSD[i] = sqrt(out$var.sig.gp.fit)
  if (i %% 1000 == 0){
    print(paste('AdaTiss fitting done ... ', 100*i/nrow(zMat),'%',sep = ''))
  }
}

write.csv(zMat,paste('additional_RNA_TS_scores_input.csv',sep = ''))
write.csv(fittingRes,paste('additional_RNA_TS_scores_populationMetrics.csv',sep = ''))

