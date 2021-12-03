# compare gap filling vs error correcting 
summaryData = readxl::read_xlsx('./output/supp_table.xlsx',sheet = 'Table S2')
table(summaryData$predicted_by_FPA)
table(summaryData$predicted_by_default_FPA)
table(summaryData$correlated)
table(summaryData$expression_type)
predType = data.frame(gapfilling = c(NA,NA), correcting = c(NA,NA),row.names = c('local','optimal'))
predType['local',"gapfilling"] = sum(summaryData$predicted_by_default_FPA=='Yes' &
  summaryData$expression_type %in% c("No GPR","Not measured") &
    summaryData$correlated == 'No')
predType['local',"correcting"] = sum(summaryData$predicted_by_default_FPA=='Yes' &
                                       !(summaryData$expression_type %in% c("No GPR","Not measured")) &
                                       summaryData$correlated == 'No')
predType["optimal","gapfilling"] = sum(summaryData$predicted_by_FPA=='Yes' &
                                       summaryData$expression_type %in% c("No GPR","Not measured") &
                                       summaryData$correlated == 'No')
predType['optimal',"correcting"] = sum(summaryData$predicted_by_FPA=='Yes' &
                                       !(summaryData$expression_type %in% c("No GPR","Not measured")) &
                                       summaryData$correlated == 'No')
# show predictions of all rxn in venn plot
library(eulerr)
dev.off()
pdf('figures/rxn_venn.pdf',width = 7,height = 7)
plot(euler(list(all = summaryData$rxnID,
                expression_only=summaryData$rxnID[summaryData$correlated == 'Yes'],
                localFPA=summaryData$rxnID[summaryData$predicted_by_default_FPA == 'Yes'],
                optimalFPA = summaryData$rxnID[summaryData$predicted_by_FPA == 'Yes'])), quantities = TRUE)
dev.off()
# classify the new predictions
dev.off()
pdf('figures/pred_type_local.pdf',width = 7,height = 7)
pie(as.numeric(predType[1,]), colnames(predType),main = paste('gapfill = ',predType[1,1], ' correcting = ',predType[1,2],sep = ''))
dev.off()
dev.off()
pdf('figures/pred_type_optimal.pdf',width = 7,height = 7)
pie(as.numeric(predType[2,]), colnames(predType),main = paste('gapfill = ',predType[2,1], ' correcting = ',predType[2,2],sep = ''))
dev.off()





PCC_all = read.csv('output/PCC_titration_all.csv',row.names = 1);
sig = read.csv('output/heatmapTbl_sigLabel.csv',row.names = 1)

boundary = read.csv('output/heatmapTbl_boundaries.csv')
colnames(PCC_all) = boundary$Var1
infoTbl = read.csv('output/summary_table_reaction_information.csv')
validRxns = infoTbl$rxnID[!infoTbl$expression_type %in% c('Not measured','No GPR')]
PCC_all = PCC_all[validRxns,]
sig_all = rownames(sig)
sig_exp = rownames(sig)[sig$sigMat2 == 1]
sig_defaultFPA = rownames(sig)[sig$sigMat1 == 1]

library(matrixStats)
PCC_all$max = rowMaxs(as.matrix(PCC_all[,3:ncol(PCC_all)]))


dev.off()
pdf('figures/PCC_comparison.pdf',width = 7,height = 7)
# plot(PCC_all$`expression only`, PCC_all$`base 2 - boundary 6`, ylim = c(-1,1), xlim = c(-1,1),
#      xlab = 'Expression only',
#      ylab = 'FPA integration',col = '#0072BD',
#      lwd = 2)
# arrows(PCC_all$`expression only`, PCC_all$`expression only`,
#        PCC_all$`expression only`, PCC_all$`base 2 - boundary 6`, col = '#0072BD',
#        length = 0.1)
# arrows(PCC_all$`expression only`, PCC_all$`base 2 - boundary 6`,
#        PCC_all$`expression only`, PCC_all$max, col = '#77AC30',
#        length = 0.1,lty = 2)
# plot for FPA
plot(PCC_all$`expression only`, PCC_all$max, col = 'black', lwd = 2,pch =16,
     xlab = 'PCC by target expression',
     ylab = 'PCC by optimal boundary integration')
# arrows(PCC_all$`expression only`, PCC_all$`expression only`,
#        PCC_all$`expression only`, PCC_all$max, col = 'black',
#        length = 0.1,lty = 2)
abline(a = 0,b = 1,lwd = 2,lty = 2)
abline(v = -0.4,lwd = 1,lty = 2)
abline(v = 0.4,lwd = 1,lty = 2)
#legend(x = -1,y = 1,legend = c('local integration','optimal boundary integration'),pch = 1,bty = 'n',lty=0,lwd = 2,col = c('#0072BD','#77AC30'), cex=1)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
box(lwd=2)
dev.off()

# distribution of the delta 
dev.off()
pdf('figures/PCC_comparison_deltaPCC.pdf',width = 7,height = 7)

# hist(PCC_all$max - PCC_all$`base 2 - boundary 6`, breaks = 30)
# hist(PCC_all$`base 2 - boundary 6` - PCC_all$`expression only`, breaks = 30)
# hist(PCC_all$max - PCC_all$`expression only`, breaks = 30)

tmp = sort(PCC_all$`expression only`)
ind = order(PCC_all$`expression only`)
plot(tmp, PCC_all$`base 2 - boundary 6`[ind] - PCC_all$`expression only`[ind], 
     ylim = c(-2,2), 
     xlim = c(-1,1),
     xlab = 'Expression only',
     ylab = 'FPA integration',col = '#0072BD',
     lwd = 2, type = 'p')
abline(a = 0,b = 0,lwd = 2,lty = 2)
# plot for FPA
tmp = sort(PCC_all$`expression only`)
ind = order(PCC_all$`expression only`)
points(tmp, PCC_all$max[ind] - PCC_all$`expression only`[ind], col = '#77AC30', lwd = 2)
legend(x = -1,y = 1,legend = c('local integration','optimal boundary integration'),pch = 1,bty = 'n',lty=0,lwd = 2,col = c('#0072BD','#77AC30'), cex=1)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
box(lwd=2)
dev.off()




