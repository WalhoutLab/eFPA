setEnvForAnalysis

%% we need to decide if we want to use common or use all for this analysis!
%% maybe should go back to comon

addpath('PlotPub/lib')
load output/FPA_rxn_protein_TS_all_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
%% compare weighted with raw
%Zcutoffs = [1.645,1.96,2.58];
z_cutoff = 1.645;
% find reactions with single gene
simpleRxns = model.rxns(sum(model.rxnGeneMat,2)==1);
% only keep rFP valid 
simpleRxns = intersect(simpleRxns,regexprep(rowlabels,'_.$','')); % if both _f and _r has rFP, _f (the previous one, is picked).

% also prepare a non-dup set 
simpleGenes = {};
for i = 1:length(simpleRxns)
    simpleGenes(i) = model.genes(model.rxnGeneMat(ismember(model.rxns,simpleRxns{i}),:)==1);
end
[simpleGenes_unique, ia] = unique(simpleGenes);
simpleRxns_unique = simpleRxns(ia);

% % load the reported Tissue-enriched and tissue-specific genes 
% enrichTbl = readtable('input/suppTbls/populationMeans.xlsx','Sheet','protein');
% gene_enrich = intersect(simpleGenes, enrichTbl.ensembl_id(ismember(enrichTbl.prt_ench_category,{'prt_enriched_not_spec','prt_specific'}))); % both enrich and specific
% gene_specific = intersect(simpleGenes, enrichTbl.ensembl_id(ismember(enrichTbl.prt_ench_category,{'prt_specific'}))); % both enrich and specific

% only look at measured genes
A = ismember(simpleGenes,TsTbl.ensembl_id);
simpleGenes = simpleGenes(A);
simpleRxns = simpleRxns(A);
A = ismember(simpleGenes_unique,TsTbl.ensembl_id);
simpleGenes_unique = simpleGenes_unique(A);
simpleRxns_unique = simpleRxns_unique(A);

% make the tissue calls 
tissueEnrichMat = zeros(length(simpleRxns),length(conditions));
for i = 1:size(tissueEnrichMat,1)
    tissueEnrichMat(i,:) = TsTbl{strcmp(TsTbl.ensembl_id,simpleGenes{i}),conditions} >= z_cutoff;
end
tissueEnrichMat_unique = zeros(length(simpleRxns_unique),length(conditions));
for i = 1:size(tissueEnrichMat_unique,1)
    tissueEnrichMat_unique(i,:) = TsTbl{strcmp(TsTbl.ensembl_id,simpleGenes_unique{i}),conditions} >= z_cutoff;
end

% the Ts Matrix
TsTbl_ori = readtable('input/suppTbls/proteinTSscore.xlsx');
[A B] = ismember(simpleGenes,TsTbl_ori.ensembl_id);
TsMat_simpleRxns = TsTbl_ori{B(A),conditions};
[A B] = ismember(simpleGenes_unique,TsTbl_ori.ensembl_id);
TsMat_simpleRxns_unique = TsTbl_ori{B(A),conditions};

%%titrate the equvalent threshold for delta rFP - by unique simple reactions 
load output/FPA_rxn_protein_TS_all_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd_local = relFP;

[A B] = ismember(simpleRxns_unique,regexprep(rowlabels,'_.$','')); % if both _f and _r has rFP, _f (the previous one, is picked).
relFP_simpleRxns_unique = relFP_wtd_local(B(A),:);
cutoffs = 0:0.01:1;
precisionMat = zeros(length(cutoffs),length(conditions));
recallMat = zeros(length(cutoffs),length(conditions));
OverallPrecision = zeros(length(cutoffs),1);
OverallRecall = zeros(length(cutoffs),1);
for i = 1: length(cutoffs)
    predictions = (relFP_simpleRxns_unique - median(relFP_simpleRxns_unique,2)) >= cutoffs(i);
    TP = sum(predictions & tissueEnrichMat_unique,1);
    FP = sum(predictions & ~tissueEnrichMat_unique,1);
    FN = sum(~predictions & tissueEnrichMat_unique,1);
    precisionMat(i,:) = TP ./ (TP + FP);
    recallMat(i,:) = TP ./ (TP + FN);
    OverallPrecision(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall(i) = sum(TP) ./ (sum(TP) + sum(FN));
end
precisionMat(isnan(precisionMat)) = 1; % 0/ 0;

% titrate the equvalent threshold for delta rFP - by unique simple reactions - raw
load output/FPA_rxn_protein_raw_all_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd_local = relFP;

[A B] = ismember(simpleRxns_unique,regexprep(rowlabels,'_.$','')); % if both _f and _r has rFP, _f (the previous one, is picked).
relFP_simpleRxns_unique = relFP_wtd_local(B(A),:);
cutoffs = 0:0.01:1;
precisionMat_raw = zeros(length(cutoffs),length(conditions));
recallMat_raw = zeros(length(cutoffs),length(conditions));
OverallPrecision_raw = zeros(length(cutoffs),1);
OverallRecall_raw = zeros(length(cutoffs),1);
for i = 1: length(cutoffs)
    predictions = (relFP_simpleRxns_unique - median(relFP_simpleRxns_unique,2)) >= cutoffs(i);
    TP = sum(predictions & tissueEnrichMat_unique,1);
    FP = sum(predictions & ~tissueEnrichMat_unique,1);
    FN = sum(~predictions & tissueEnrichMat_unique,1);
    precisionMat_raw(i,:) = TP ./ (TP + FP);
    recallMat_raw(i,:) = TP ./ (TP + FN);
    OverallPrecision_raw(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall_raw(i) = sum(TP) ./ (sum(TP) + sum(FN));
end
precisionMat_raw(isnan(precisionMat_raw)) = 1; % 0/ 0;

figure(5)
hold on
plot(OverallRecall, OverallPrecision,'-');
plot(OverallRecall_raw, OverallPrecision_raw,'-');
plot(OverallRecall, cutoffs,'-');
plot(OverallRecall_raw, cutoffs,'-');
hold off
legend({'PR-curve (shrunk)','PR-curve (raw)','cutoff-recall (shrunk)','cutoff-recall (raw)'});
ylabel('Precison (or cutoff)')
xlabel('recall')
xlim([0 1]);
ylim([0 1]);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [6, 6];
plt.LineWidth = 3;
plt.FontSize = 15;
plt.XTick = plt.YTick;
plt.FontName = 'Arial';
% plt.export('figures/benchmarkTSweighting.svg');
%saveas(figure(5), 'figures/benchmarkTSweighting.svg');

figure(6)
hold on
F1 = (2 .* OverallPrecision .* OverallRecall) ./ (OverallPrecision + OverallRecall);
F2 = (2 .* OverallPrecision_raw .* OverallRecall_raw) ./ (OverallPrecision_raw + OverallRecall_raw);
plot(cutoffs,F1,'.-');
plot(cutoffs,F2,'.-');
hold off
legend({'shrunk','raw'});
xlabel('cutoff')
ylabel('F measure')
%% correlation with z-score = not very useful and hard to explain 
%% titrate the equvalent threshold for delta rFP - by unique simple reactions 
load output/FPA_rxn_protein_TS_all_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd_local = relFP;
[A B] = ismember(simpleRxns_unique,regexprep(rowlabels,'_.$','')); % if both _f and _r has rFP, _f (the previous one, is picked).
relFP_simpleRxns_unique = relFP_wtd_local(B(A),:);
tmp = -log10(2.*(1-normcdf(abs(TsMat_simpleRxns_unique')))) .* sign(TsMat_simpleRxns_unique');
PearsonCorr = diag(corr(tmp,relFP_simpleRxns_unique','type','Pearson'));

load output/FPA_rxn_protein_raw_all_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd_local = relFP;

[A B] = ismember(simpleRxns_unique,regexprep(rowlabels,'_.$','')); % if both _f and _r has rFP, _f (the previous one, is picked).
relFP_simpleRxns_unique = relFP_wtd_local(B(A),:);
tmp = -log10(2.*(1-normcdf(abs(TsMat_simpleRxns_unique')))) .* sign(TsMat_simpleRxns_unique');
PearsonCorr_raw = diag(corr(tmp,relFP_simpleRxns_unique','type','Pearson'));
figure
[A B] = sort(PearsonCorr_raw,'descend');
hold on 
bar(A);
bar(PearsonCorr(B));
hold off
legend({'raw','shrunk'});




