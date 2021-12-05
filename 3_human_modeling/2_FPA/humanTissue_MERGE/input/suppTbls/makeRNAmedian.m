allSamples = readtable('rnaAllSamples_log2TPM.xlsx');
sampleSheet = readtable('rnaAllSamples_log2TPM.xlsx','Sheet','sampleSheet');
% remove duplicates (unclear what they are)
allSamples = allSamples(cellfun(@(x) isempty(regexp(x,'_PAR_Y$','once')),allSamples.gene_id_full),:);
%% only look at those with TS scores
% TsTbl = readtable('RNATSscore.xlsx');
% allSamples = allSamples(ismember(allSamples.gene_id,TsTbl.ensembl_id),:);
%%
allTissues = unique(sampleSheet.tissueID);
tissueMedians = nan(length(allSamples.gene_id_full),length(allTissues));
for i = 1:length(allTissues)
    sampleIDs = sampleSheet.RNA_sample_ID(strcmp(sampleSheet.tissueID,allTissues{i}));
    sampleIDs = regexprep(sampleIDs,'-','_');
    tissueMedians(:,i) = median(allSamples{:,sampleIDs},2);
end
tissueMedians = array2table(tissueMedians);
tissueMedians.Properties.RowNames = allSamples.gene_id;
tissueMedians.Properties.VariableNames = regexprep(allTissues,' ','');
writetable(tissueMedians,'RNATissueMedian_log2TPM_homemade.xlsx','FileType','spreadsheet','WriteRowNames',true,'WriteVariableNames',true);
%% make the corresponding RNA-TS score 
% we still use the original TS score if avaiable 
allSamples = readtable('additional_RNA_TS_scores_input.csv');
sampleSheet = readtable('rnaAllSamples_log2TPM.xlsx','Sheet','sampleSheet');
% remove duplicates (unclear what they are)
allSamples = allSamples(cellfun(@(x) isempty(regexp(x,'_PAR_Y$','once')),allSamples.Var1),:);
%% 
allTissues = unique(sampleSheet.tissueID);
tissueMedians = nan(length(allSamples.Var1),length(allTissues));
for i = 1:length(allTissues)
    sampleIDs = sampleSheet.RNA_sample_ID(strcmp(sampleSheet.tissueID,allTissues{i}));
    sampleIDs = regexprep(sampleIDs,'-','_');
    tissueMedians(:,i) = median(allSamples{:,sampleIDs},2);
end
tissueMedians = array2table(tissueMedians);
IDtbl = readtable('additional_RNA_TS_IDtbl.csv');
[A B] = ismember(allSamples.Var1,IDtbl.data_gene_id_full);
tissueMedians.Properties.RowNames = IDtbl.data_gene_id(B(A));
tissueMedians.Properties.VariableNames = regexprep(allTissues,' ','');
writetable(tissueMedians,'RNATS_score_homemade.xlsx','FileType','spreadsheet','WriteRowNames',true,'WriteVariableNames',true);

fitting = readtable('additional_RNA_TS_scores_populationMetrics.csv');
fitting = fitting(ismember(fitting.geneID, allSamples.Var1),:);
[A B] = ismember(fitting.geneID,IDtbl.data_gene_id_full);
fitting.geneID = IDtbl.data_gene_id(B(A));
writetable(fitting,'RNATS_fitting_homemade.xlsx','FileType','spreadsheet','WriteRowNames',true,'WriteVariableNames',true);

%% compare the homemade TS with the published 
pubTS = readtable('RNATSscore.xlsx');
pubTS.Properties.VariableNames(14) = {'GEjunction'};
pubTS.Properties.RowNames = pubTS.ensembl_id;
pubTS = pubTS(:,5:end);
homeTS = tissueMedians(pubTS.Properties.RowNames, pubTS.Properties.VariableNames);
%%
plot([homeTS{:,:}], [pubTS{:,:}],'.')
