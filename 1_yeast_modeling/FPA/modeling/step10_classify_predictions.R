# About: analyze the type of prediction and benefit from integration

# compare gap filling vs error correcting 
summaryData = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
table(summaryData$predicted_by_optimal_boundary_FPA)
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
predType["optimal","gapfilling"] = sum(summaryData$predicted_by_optimal_boundary_FPA=='Yes' &
                                       summaryData$expression_type %in% c("No GPR","Not measured") &
                                       summaryData$correlated == 'No')
predType['optimal',"correcting"] = sum(summaryData$predicted_by_optimal_boundary_FPA=='Yes' &
                                       !(summaryData$expression_type %in% c("No GPR","Not measured")) &
                                       summaryData$correlated == 'No')
# show predictions of all rxn in venn plot
library(eulerr)
dev.off()
pdf('figures/rxn_venn.pdf',width = 7,height = 7)
plot(euler(list(all = rownames(summaryData),
                expression_only=rownames(summaryData)[summaryData$correlated == 'Yes'],
                localFPA=rownames(summaryData)[summaryData$predicted_by_default_FPA == 'Yes'],
                optimalFPA = rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA == 'Yes'])), quantities = TRUE)
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

# directly compare the PCC for understanding the benefit
PCC_all = read.csv('output/PCC_titration_all.csv',row.names = 1);
boundary = data.frame(Var1 = c(c('base2-boundary6','expression only'), seq(0,40,0.5)))
colnames(PCC_all) = boundary$Var1
infoTbl = read.csv('output/summary_table_reaction_information.csv')
validRxns = infoTbl$rxnID[!infoTbl$expression_type %in% c('Not measured','No GPR')]
PCC_all = PCC_all[validRxns,]
library(matrixStats)
PCC_all$max = rowMaxs(as.matrix(PCC_all[,3:ncol(PCC_all)]))
#plot
dev.off()
pdf('figures/PCC_comparison.pdf',width = 7,height = 7)
plot(PCC_all$`expression only`, PCC_all$max, col = 'black', lwd = 2,pch =16,
     xlab = 'PCC by target expression',
     ylab = 'PCC by optimal boundary integration')
abline(a = 0,b = 1,lwd = 2,lty = 2)
abline(v = -0.4,lwd = 1,lty = 2)
abline(v = 0.4,lwd = 1,lty = 2)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
box(lwd=2)
dev.off()





