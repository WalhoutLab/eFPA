library(pheatmap)
library(RColorBrewer)
# the annotation of reactions
# the manual connected pathway was used; pathways with reactions less than 3 were grouped to others; isolated transporters are labeled out (inner and exchange)
rMat = read.csv('output/crossCorrelation_rMat.csv',row.names = 1)
FDRmat = read.csv('output/crossCorrelation_FDRMat.csv',row.names = 1)

heatTbl = rMat
# row-wise normalize 
heatTbl = t(scale(t(heatTbl),scale = apply(heatTbl, 1, max),center = F))
sigTbl = ifelse(rMat > 0 & FDRmat < 0.05, '+','')
dev.off()
pdf('figures/NoTrack_cross_correlation.pdf',width = 60,height = 60)
# add annotation
library(xlsx)
annTbl = read.xlsx('./pathway_annotations.xlsx','more_info',header = T)
rownames(annTbl) = annTbl$rxn
annotationList = c('Purine metabolism','Phenylalanine, tyrosine and tryptophan biosynthesis','Histidine metabolism','Lysine biosynthesis',
                   'Arginine biosynthesis','Threonine, methionine and cysteine synthesis','Proline biosynthesis','Pyrimidine metabolism',
                   'Glycolysis','TCA cycle', 'Mannan synthesis',
                   'UDP-D-glucose metabolism','Fatty acid biosynthesis',
                   'transporter [inner]','transporter [exchange]','Others')
annTbl$myAnnotation = annTbl$manual_pathway
annTbl$myAnnotation[!(annTbl$myAnnotation %in% annotationList)] = 'Others'
colors = read.csv('figures/simpsons_color.csv',header = F)
colorList2 = colors$V1[2:(1+length(annotationList))]
names(colorList2) = annotationList
colorList2 = list(myAnnotation = colorList2)
colorList2$myAnnotation = colorList2$myAnnotation[unique(annTbl[,c("myAnnotation")])]
# assign annotations
ann2 = annTbl[rownames(heatTbl),"manual_pathway"]
labels_row0 = paste(rownames(heatTbl),ann2)
ann3 = annTbl[colnames(heatTbl),"manual_pathway"]
labels_col0 = paste(colnames(heatTbl),ann3)
# also add the PCC values (max PCC for each row) in the annotation
# PCCs = read.csv('output/PCC_titration_all.csv',row.names = 1)
# maxPCC = rowMaxs(as.matrix(PCCs[rownames(annTbl),3:ncol(PCCs)]))
# labels_row0 = paste(labels_row0,round(maxPCC,2))

pheatmap(heatTbl, breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                                      "RdYlBu")))(1000),
         annotation_row = annTbl[,c("myAnnotation"),drop = F],
         annotation_colors = colorList2, labels_row = labels_row0, labels_col = labels_col0,
         border_color = 0,
         fontsize = 6,fontsize_row = 6, fontsize_col = 6,
         cellwidth = 6, cellheight = 6,
         cluster_rows = T,cluster_cols = T,clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation',
         display_numbers = sigTbl, number_color = 'Black',fontsize_number = 3
)

dev.off()


# compare the best PCC here with the best PCC in optimal bound FPA
library(matrixStats)
boundary = data.frame(Var1 = c(c('base2-boundary6','expression only'), seq(0,40,0.5)))
FPAtbl = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
PCCs = read.csv('output/PCC_titration_all_simpleDecay.csv',row.names = 1)
colnames(PCCs) = boundary$Var1
FPAtbl$maxPCC_FPA = rowMaxs(as.matrix(PCCs[rownames(FPAtbl),3:ncol(PCCs)]))
FPAtbl$maxPCC_all2all = rowMaxs(as.matrix(rMat[rownames(FPAtbl),]))
FPAtbl$delta = FPAtbl$maxPCC_FPA - FPAtbl$maxPCC_all2all

plot(FPAtbl$maxPCC_FPA, FPAtbl$maxPCC_all2all, xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b= 1)
# ==> FPA underperform for most rxns, except for some transporter that is related to the media settings


# plot the heatmap only for correlated rxns against others
rMat = read.csv('output/coCorrelation_rMat.csv',row.names = 1)
heatTbl = rMat
# row-wise normalize: row is each ROI flux under consideration; col is each ROI expression
maxCorr = rowMaxs(as.matrix(heatTbl))
# clean up
#heatTbl = heatTbl[,!colAlls(heatTbl==0)]
#heatTbl = heatTbl[maxCorr > 0.12,] # 0.4*0.4*0.75
#maxCorr = maxCorr[maxCorr > 0.12]
heatTbl = t(scale(t(heatTbl),scale = apply(heatTbl, 1, max),center = F))
#heatTbl = heatTbl[,colMaxs(heatTbl) > 0.5]
#write.csv(rownames(heatTbl),file = 'output/indicatorPoints.csv')

# the rows will be ordered by flux
fluxMat = read.csv('output/supp1B_normalizedFlux.csv',row.names = 1)
fluxMat = fluxMat[rownames(rMat),]
rowOrder = as.dist(1 - cor(t(fluxMat), use = "pairwise.complete.obs", method = "pearson"))

dev.off()
pdf('figures/NoTrack_co_correlation_stringent.pdf',width = 60,height = 60)
# add annotation
library(xlsx)
library(grid)
annTbl = read.xlsx('./pathway_annotations.xlsx','more_info',header = T)
summaryData = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
annTbl$predicted_FPA = ifelse(annTbl$rxn %in% rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA ==
                                                                      'Yes'], 'Yes','No')
rownames(annTbl) = annTbl$rxn
colorList2 = c('#77AC30','#000000')
names(colorList2) = c('Yes','No')
colorList2 = list(predicted_FPA = colorList2)
# assign annotations
ann2 = annTbl[rownames(heatTbl),"manual_pathway"]
labels_row0 = paste(rownames(heatTbl),ann2)
ann3 = annTbl[colnames(heatTbl),"manual_pathway"]
labels_col0 = paste(colnames(heatTbl),ann3)
# also add the PCC values (max PCC for each row) in the annotation
labels_row0 = paste(labels_row0,round(maxCorr,2))

heatTbl = as.data.frame(heatTbl)
heatTbl$colors=ifelse(rownames(heatTbl) %in% colnames(heatTbl),"red","black")
e = pheatmap(heatTbl[,1:(ncol(heatTbl)-1)], breaks = seq(0,1,0.001),color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                                    "RdYlBu")))(1000),
         annotation_row = annTbl[,c("predicted_FPA"),drop = F],
         annotation_colors = colorList2, labels_row = labels_row0, labels_col = labels_col0,
         border_color = 0,
         fontsize = 6,fontsize_row = 6, fontsize_col = 6,
         cellwidth = 6, cellheight = 6,
         cluster_rows = T,cluster_cols = T
)
cols=heatTbl[order(match(labels_row0, e$gtable$grobs[[5]]$label)), ]$colors
e$gtable$grobs[[5]]$gp=gpar(col=cols)
print(e)
dev.off()

# new fluxes
allFluxes = union(rownames(heatTbl),rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA =='Yes'])
length(allFluxes)/232
length(intersect(rownames(heatTbl),rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA =='Yes']))


# make the new venn diagram
summaryData = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
# show predictions of all rxn in venn plot
library(eulerr)
dev.off()
pdf('figures/rxn_venn_with_indicatorRxns.pdf',width = 7,height = 7)
plot(euler(list(all = rownames(summaryData),
                expression_only=rownames(summaryData)[summaryData$correlated == 'Yes'],
                localFPA=rownames(summaryData)[summaryData$predicted_by_default_FPA == 'Yes'],
                optimalFPA = rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA == 'Yes'],
                indicatorRxns = rownames(heatTbl))), quantities = TRUE)
dev.off()



# lastly (maybe?), plot the pathway annotation of diff type of dictation
rMat = read.csv('output/coCorrelation_rMat.csv',row.names = 1)
dictated_rxns = rownames(rMat)
summaryData = read.csv('output/summary_table_reaction_information.csv',row.names = 1)
pathwayInfo = xlsx::read.xlsx('pathway_annotations.xlsx','more_info')
allCount = table(pathwayInfo$manual_pathway)
allCount = as.data.frame(allCount)
# correlated 
set1 = rownames(summaryData)[summaryData$correlated == 'Yes']
correlatedRxns = pathwayInfo$manual_pathway[pathwayInfo$rxn %in% set1]
correlatedRxnsCount = table(correlatedRxns)
allCount$correlated = correlatedRxnsCount[as.character(allCount$Var1)]
# pathway-level predicted rxns
set2 = rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA == 'Yes']
correlatedRxns = pathwayInfo$manual_pathway[pathwayInfo$rxn %in% setdiff(set2,set1)]
correlatedRxnsCount = table(correlatedRxns)
allCount$FPApredicted = correlatedRxnsCount[as.character(allCount$Var1)]
# indicator predicted rxns
correlatedRxns = pathwayInfo$manual_pathway[pathwayInfo$rxn %in% setdiff(dictated_rxns,union(set1,set2))]
correlatedRxnsCount = table(correlatedRxns)
allCount$dictated = correlatedRxnsCount[as.character(allCount$Var1)]
# not predicted
myrxn = setdiff(rownames(summaryData),union(dictated_rxns, rownames(summaryData)[summaryData$predicted_by_optimal_boundary_FPA == 'Yes']))
correlatedRxns = pathwayInfo$manual_pathway[pathwayInfo$rxn %in% myrxn]
correlatedRxnsCount = table(correlatedRxns)
allCount$not_predicted = correlatedRxnsCount[as.character(allCount$Var1)]
# clean up
allCount[is.na(allCount)] = 0
allCount$Var1 = factor(as.character(allCount$Var1),levels = allCount$Var1[order(allCount$Freq)])
countMat = allCount[,c(3,4,5,6)]
rownames(countMat) = allCount$Var1
countMat = reshape2::melt(as.matrix(countMat))
countMat$Var1 = factor(as.character(countMat$Var1),levels = allCount$Var1[order(allCount$Freq)])
countMat$Var2 = factor(as.character(countMat$Var2),levels = c('not_predicted','dictated','FPApredicted','correlated'))
library(artyfarty)
library(ggplot2)
p = ggplot(countMat, aes(x= Var1, y=value, fill = Var2))+
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
        panel.border = element_blank())+
  scale_fill_manual("legend", values = c("not_predicted" = "#C8C2BC", 
                                         "correlated" = "#FF6701", 
                                         "FPApredicted" = "#FEA82F",
                                         "dictated" = "#FFC288"))

dev.off()
pdf('figures/pathway_annotation_all_rxns.pdf',width = 5,height = 3.25)
print(p)
dev.off()

