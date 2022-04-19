fluxMat = read.csv('output/flux_expression_correlation_flux_mat_to_plot.csv')
expMat = read.csv('output/flux_expression_correlation_exp_mat_to_plot.csv')
sigInfo = read.csv('output/summary_table_reaction_information.csv')
pdf('figures/correlationPlots/merged.pdf',width = 8.5,height = 11)
par(mfrow=c(13,12),mar=c(0,0.1,1,0.1), oma = c(3,3,5,0))
for (i in 1:156){
  pcc = cor(expMat[,i], fluxMat[,i])
  txt = paste(colnames(expMat)[i],', PCC=', round(pcc,2),sep = '')
  col = ifelse(colnames(expMat)[i] %in% sigInfo$rxnID[sigInfo$correlated== 'Yes'],
               'red','black')
  xlim = c(min(expMat[,i]) - 0.1, max(expMat[,i]) + 0.1)
  ylim = c(min(fluxMat[,i]) - (0.1*max(fluxMat[,i])), max(fluxMat[,i])*1.1)
  
  plot(expMat[,i], fluxMat[,i],pch = 20,cex = 0.6,
       ylab="",yaxt="n", xlab="",xaxt="n",
       xlim = xlim, ylim = ylim
       )
  #lines(expMat[,i], predict(lm(fluxMat[,i]~expMat[,i])),col='grey', lty = 1)
  box(col=col)
  title(txt, line = 0.2, cex.main = 0.6)
}
dev.off()



