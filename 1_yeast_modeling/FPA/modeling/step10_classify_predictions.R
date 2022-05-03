# About: visualization of correlated rxns, analyze the type of prediction and benefit from integration

# pathway annotation of correlated rxns
# also label for each pathway how many are tested
summaryData = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
pathwayInfo = xlsx::read.xlsx('pathway_annotations.xlsx','more_info')
correlatedRxns = pathwayInfo$manual_pathway[pathwayInfo$rxn %in% rownames(summaryData)[summaryData$correlated == 'Yes']]
correlatedRxnsCount = table(correlatedRxns)
correlatedRxnsCount = as.data.frame(correlatedRxnsCount)
allCount = table(pathwayInfo$manual_pathway[pathwayInfo$rxn %in% rownames(summaryData)[summaryData$expression_type %in% c('Low Diversity','High Diversity','Low Variability')]])
correlatedRxnsCount$allCount = allCount[as.character(correlatedRxnsCount$correlatedRxns)]
correlatedRxnsCount$correlatedRxns = factor(as.character(correlatedRxnsCount$correlatedRxns),levels = correlatedRxnsCount$correlatedRxns[order(correlatedRxnsCount$allCount)])
correlatedRxnsCount$not_correlated = correlatedRxnsCount$allCount-correlatedRxnsCount$Freq
correlatedRxnsCount = correlatedRxnsCount[,c(1,2,4)]
colnames(correlatedRxnsCount) = c('pathway','correlated','not correlated')
correlatedRxnsCount = reshape2::melt(correlatedRxnsCount)
correlatedRxnsCount$variable = factor(as.character(correlatedRxnsCount$variable),levels = c('not correlated','correlated'))
# add the annotation of enrichment for each pathway 
# assessible reaction
validRxns = rownames(summaryData)[summaryData$expression_type %in% c('Low Diversity','High Diversity','Low Variability')]
allPath = unique(pathwayInfo$manual_pathway[pathwayInfo$rxn %in% validRxns])
allHits = sum(summaryData$correlated == 'Yes')
pvalues = c()
for (i in 1:length(allPath)){
  hits = length(intersect(pathwayInfo$rxn[pathwayInfo$manual_pathway == allPath[i]], 
                          rownames(summaryData)[summaryData$correlated == 'Yes']))
  allIn = length(intersect(pathwayInfo$rxn[pathwayInfo$manual_pathway == allPath[i]], 
                           validRxns))
  pvalues[i] = phyper(hits-1, allIn, length(validRxns) - allIn,allHits, lower.tail = F)

}
pvalues = p.adjust(pvalues,method = 'BH')
names(pvalues) = allPath
correlatedRxnsCount$padj = pvalues[as.character(correlatedRxnsCount$pathway)]

library(artyfarty)
library(ggplot2)
p = ggplot(correlatedRxnsCount, aes(x= pathway, y=value, fill = variable))+
  geom_bar(stat="identity",position="stack", size=0.25,colour="black")+
  ylab("Number of reactions")+
  coord_flip()+ scale_y_continuous(expand = c(0, 0)) +
  theme_bw()+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(colour ='black'),axis.text.y = element_text(colour ='black'),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #legend.title = element_blank(),
        axis.line.y.left = element_line(colour = 'black', size = 0.5),
        axis.line.x.bottom = element_line(colour = 'black', size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_blank())
dev.off()
pdf('figures/pathway_annotation_correlated_rxns.pdf',width = 5,height = 2)
print(p)
dev.off()
# same for other predictions
correlatedRxns = pathwayInfo$manual_pathway[pathwayInfo$rxn %in% rownames(summaryData)[summaryData$predicted_by_default_FPA == 'Yes']]
correlatedRxnsCount = table(correlatedRxns)
correlatedRxnsCount = as.data.frame(correlatedRxnsCount)
library(artyfarty)
library(ggplot2)
p = ggplot(correlatedRxnsCount, aes(x=reorder(correlatedRxns, -Freq), y=Freq))+
  geom_bar(stat="identity",position="identity", size=0.25,colour="black")+
  ylab("Number of reactions")+
  coord_flip()+ scale_y_continuous(expand = c(0, 0)) +
  theme_bw()+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(colour ='black'),axis.text.y = element_text(colour ='black'),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.line.y.left = element_line(colour = 'black', size = 0.5),
        axis.line.x.bottom = element_line(colour = 'black', size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_blank())
dev.off()
pdf('figures/pathway_annotation_predicted_by_defaultFPA_rxns.pdf',width = 4,height = 2.75)
print(p)
dev.off()

correlatedRxns = pathwayInfo$manual_pathway[pathwayInfo$rxn %in% rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA == 'Yes']]
correlatedRxnsCount = table(correlatedRxns)
correlatedRxnsCount = as.data.frame(correlatedRxnsCount)
library(artyfarty)
library(ggplot2)
p = ggplot(correlatedRxnsCount, aes(x=reorder(correlatedRxns, -Freq), y=Freq))+
  geom_bar(stat="identity",position="identity", size=0.25,colour="black")+
  ylab("Number of reactions")+
  coord_flip()+ scale_y_continuous(expand = c(0, 0)) +
  theme_bw()+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(colour ='black'),axis.text.y = element_text(colour ='black'),
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.line.y.left = element_line(colour = 'black', size = 0.5),
        axis.line.x.bottom = element_line(colour = 'black', size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_blank())
dev.off()
pdf('figures/pathway_annotation_predicted_by_optimalBondFPA_rxns.pdf',width = 4,height = 3.25)
print(p)
dev.off()

# compare gap filling vs error correcting 
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





