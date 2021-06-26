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