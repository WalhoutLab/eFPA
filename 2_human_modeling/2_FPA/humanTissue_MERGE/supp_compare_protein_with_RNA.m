% load data 
addpath ../scripts/

% read the FPA results
setEnvForAnalysis

% compare the TS score 
TsTbl_RNA = readtable('input/suppTbls/RNATSscore.xlsx');
TsTbl_pro = readtable('input/suppTbls/proteinTSscore.xlsx');
commGenes = intersect(TsTbl_RNA.ensembl_id, TsTbl_pro.ensembl_id);
TsTbl_pro.Properties.RowNames = TsTbl_pro.ensembl_id;
TsTbl_RNA.Properties.RowNames = TsTbl_RNA.ensembl_id;

Tbl1 = table2array(TsTbl_pro(commGenes, 5:36));
Tbl2 = table2array(TsTbl_RNA(commGenes, 5:36));

r = [];
for i = 1:size(Tbl1)
   r(i) = corr(Tbl1(i,:)', Tbl2(i,:)','Rows','complete');
end
figure
h1 = histogram(r, 'Normalization','probability');

% plot metabolic genes 
hold on
histogram(r(ismember(commGenes, model.genes)), 'Normalization','probability', 'BinEdges', h1.BinEdges);

% plot the rFP correlation 

% load FPA protein and RNA res dataset
load output/FPA_rxn_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;

load output/FPA_rxn_RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
relFP(rmInd,:) = [];
relFP_RNA = relFP;

% remove numeric issues
pass = ~any(abs(relFP_wtd)>1.01,2) & ~any(abs(relFP_RNA)>1.01,2);
rowlabels = rowlabels(pass);
relFP_RNA = relFP_RNA(pass,:);
relFP_wtd = relFP_wtd(pass,:);

% we dont filter for visualization (cutoff = 0)
relFP_filtered = relFP_wtd - median(relFP_wtd,2);
rowlabels_filtered = rowlabels;
relFP_filtered_RNA = relFP_RNA - median(relFP_RNA,2);
keep = any(abs(relFP_filtered) > 0 | abs(relFP_filtered_RNA) > 0,2);
relFP_filtered = relFP_filtered(keep,:);
rowlabels_filtered = rowlabels_filtered(keep);
relFP_filtered_RNA = relFP_filtered_RNA(keep,:);


r2 = [];
for i = 1:size(relFP_filtered)
   r2(i) = corr(relFP_filtered(i,:)', relFP_filtered_RNA(i,:)','Rows','complete');
end

histogram(r2, 'Normalization','probability', 'BinEdges', h1.BinEdges);
legend({'relative expression - all genes','relative expression - metabolic genes','rFP - all reactions'})
hold off

xlabel('Pearson correlation coefficient')
ylabel('Probability')


addpath('PlotPub/lib')
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [4, 3.5];
plt.LineWidth = 1;
plt.FontSize = 12;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Interpreter = 'none';
plt.TickDir = 'out';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.export(['figures/comp_RNA_vs_protein_data.pdf']);



% double check using the raw expression (logTPM and logFC)
% FcTbl = readtable('input/suppTbls/proteinTissueMedian_log2FC_againts_reference.xlsx');
% logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
% commGenes = intersect(FcTbl.gene_id, logTPMTbl.gene_id);
% FcTbl.Properties.RowNames = FcTbl.gene_id;
% logTPMTbl.Properties.RowNames = logTPMTbl.gene_id;
% 
% Tbl1 = table2array(FcTbl(commGenes, TsTbl_pro.Properties.VariableNames(5:36)));
% Tbl2 = table2array(TsTbl_RNA(commGenes, TsTbl_pro.Properties.VariableNames(5:36)));
% 
% r = [];
% for i = 1:size(Tbl1)
%    r(i) = corr(Tbl1(i,:)', Tbl2(i,:)','Rows','complete');
% end
% figure
% histogram(r)
% results are similar - just go with TS score because it is Z-score scale
% and easy to explain 