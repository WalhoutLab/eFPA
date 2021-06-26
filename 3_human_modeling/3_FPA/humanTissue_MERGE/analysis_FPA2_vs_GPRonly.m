setEnvForAnalysis
%%
load output/FPA_rxn_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
          
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;

load output/FPA_rxn_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
relFP(rmInd,:) = [];
relFP_GPR = relFP;
%% thresholding
% cutoffs = 0:0.01:1;
% for i = 1: length(cutoffs)
%     N_pass_rna(i) = sum(sum(abs(relFP_GPR - median(relFP_GPR,2)) > cutoffs(i)));
%     N_pass_wtd(i) = sum(sum(abs(relFP_wtd - median(relFP_wtd,2)) > cutoffs(i)));
% end
% figure()
% hold on
% plot(cutoffs,(N_pass_rna))
% plot(cutoffs,(N_pass_wtd))
% legend({'log2(rna+1)','log2(protein+1)'});
% hold off
% xlabel('delta rFP cutoff');
% ylabel('number of tissue enrichment calls');
% %% check the most different rxns
% relFP_raw_centered = relFP_raw - median(relFP_raw,2);
% relFP_wtd_centered = relFP_wtd - median(relFP_wtd,2);
% [A B] = sort(sum(abs(relFP_raw_centered - relFP_wtd_centered),2),'descend');
% highestDiff = rowlabels(B);
% highestDiff(1:50)

%% filtering == use median and find a good cutoff
%% FC cutoff
relFP_filtered = relFP_wtd - median(relFP_wtd,2);
rowlabels_filtered = rowlabels;
relFP_filtered_GPR = relFP_GPR - median(relFP_GPR,2);

keep = any(abs(relFP_filtered) > 0.2 | abs(relFP_filtered_GPR) > 0.2,2);

relFP_filtered = relFP_filtered(keep,:);
rowlabels_filtered = rowlabels_filtered(keep);
relFP_filtered_GPR = relFP_filtered_GPR(keep,:);

length(rowlabels_filtered)
%% check for OFD consistency
lib = 'protein';
for i = 1:length(conditions)
    load(['./../../1_iMAT++/output/humanModel/',lib,'/OFD/',conditions{i},'.mat']);
    % load(['./../1_iMAT++/output/NHR_005side_speedMode2_greedyCat/',lib,'/FVA/',conditions{i},'_levels.mat']);
    eval([conditions{i},'=myCSM;']);
    % eval([conditions{i},'_levels_f=levels_f;']);
    % eval([conditions{i},'_levels_r=levels_r;']);
end
% select by high sites in OFD
IsInOFD = false(length(rowlabels_filtered),1);
for j = 1:length(rowlabels_filtered)
    highestSites = conditions((relFP_filtered(j,:) > 0.2));
    for i = 1:length(highestSites)
        eval(['myCSM=',highestSites{i},';']);
    %     eval(['levels_f=',highestStrain,'_levels_f;']);
    %     eval(['levels_r=',highestStrain,'_levels_r;']);
        myrxn = rowlabels_filtered{j};
        myrxnID = myrxn(1:end-2);
        if any(strcmp(myrxnID,targetRxns)) % reaction centric
            if strcmp(myrxn(end),'f')
                if myCSM.OFD(strcmp(model.rxns,myrxnID)) > 0
                    IsInOFD(j) = true;
                    break
                end
            else
                if myCSM.OFD(strcmp(model.rxns,myrxnID)) < 0
                    IsInOFD(j) = true;
                    break
                end
            end
        else % metabolite centric
            if any(strcmp(model.rxns,myrxnID)) % existed rxns
                if strcmp(myrxn(end),'f')
                    if levels_f(strcmp(model.rxns,myrxnID)) >= 0
                        IsInOFD(j) = true;
                        break
                    end
                else
                    if levels_r(strcmp(model.rxns,myrxnID)) >= 0
                        IsInOFD(j) = true;
                        break
                    end
                end
            else % newly added metabolite demand; not comparable with OFD
                IsInOFD(j) = true;
                break
            end
        end
    end
end
IsInOFD_pro = IsInOFD;

sum(IsInOFD_pro)
% we skip the comparisons on reactions that are not in OFD. Most not OFD
% reactions are one direction of a reversibe reaction; some are
% network-based prediction that low reactions are pushed to zero. They are
% hard to distinguish automatically, so we avoid comparision on these
% reactions
relFP_filtered = relFP_filtered(IsInOFD_pro,:);
relFP_filtered_GPR = relFP_filtered_GPR(IsInOFD_pro,:);
rowlabels_filtered = rowlabels_filtered(IsInOFD_pro,:);

%% plot relative FP
distMethod = 'cosine';
IDs = regexprep(conditions,'_','-');
IDs = [cellfun(@(x) [x,'-pro'],IDs,'UniformOutput',0),cellfun(@(x) [x,'-GPR'],IDs,'UniformOutput',0)];

relFP_filtered2 = [relFP_filtered,relFP_filtered_GPR];%normalize(relFP_filtered,2,'scale','mad');

cgo=clustergram(relFP_filtered2,'RowLabels',rowlabels_filtered,'ColumnLabels',IDs,'RowPDist',distMethod,'ColumnPDist',distMethod);
c=get(cgo,'ColorMap');
n = 100;
tmp = [ones(n,1), linspace(1,0,n)',linspace(1,0,n)'];
tmp = tmp(2:end,:);
cpr=[linspace(0,1,n)',linspace(0,1,n)',ones(n,1);...
    tmp];
set(cgo,'ColorMap',cpr);
set(cgo,'Symmetric',true);

set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
%% save figure
S = hgexport('readstyle','default_sci');
style.Format = 'svg';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
set(gca,'FontSize',15);
%%
set(h, 'FontSize', 10)
saveas(gca,'figures/clustergram_rxn_vsGPRonly_FPA.tiff');

%%
ROI = mygroup.RowNodeNames;
ROI = regexprep(ROI,'_.$','');
[ROI_annotated, ROI_enrichment] = annotateRxnSet(ROI,model);
%% also check a rxn
rxn = 'HMR_3477_f';
figure(4)
c = categorical(conditions);
bar(c, relFP_raw(strcmp(rowlabels,rxn),:))
ylabel('raw')
ylim([0,1])
figure(5)
c = categorical(conditions);
bar(c, relFP_wtd(strcmp(rowlabels,rxn),:))
ylabel('weighted')
ylim([0,1])
%% also check a rxn
rxn = 'HMR_3477_f';
figure(6)
c = categorical(conditions);
bar(c, normalize(relFP_raw(strcmp(rowlabels,rxn),:),'center','median'))
ylabel('raw')
ylim([-0.5,1])
figure(7)
c = categorical(conditions);
bar(c, normalize(relFP_wtd(strcmp(rowlabels,rxn),:),'center','median'))
ylabel('weighted')
ylim([-0.5,1])
printGPRForRxns(model,rxn(1:end-2));
%% check a gene expression
gene = 'SDS';
name2 = TsTbl.ensembl_id{strcmp(TsTbl.hgnc_symbol,gene)};
%%
name2 = 'ENSG00000142920';
%expTbl(strcmp(expTbl.ensembl_id,name2),:)
figure(1)
c = categorical(conditions);
bar(c, expTbl_shrunk{strcmp(expTbl_shrunk.ensembl_id,name2),5:end})
ylabel('fold change (shrunk, against population mean)')

figure(2)
c = categorical(conditions);
bar(c, 2.^(FcTbl{strcmp(FcTbl.gene_id,name2),conditions}))
ylabel('fold change (raw, against meta-reference)')

figure(3)
c = categorical(conditions);
bar(c, TsTbl{strcmp(TsTbl.ensembl_id,name2),conditions})
ylabel('TS scores')

figure(8)
c = categorical(conditions);
bar(c, expTbl_raw{strcmp(expTbl_raw.ensembl_id,name2),conditions})
ylabel('fold change (raw, against population mean)')

%%
ROI = mygroup.RowNodeNames;
ROI = regexprep(ROI,'_.$','');
[ROI_annotated, ROI_enrichment] = annotateRxnSet(ROI,model);

%% titrate the proper cutoff
figure(3)
figure(4)
Zcutoffs = [1.645,1.96,2.58];
for zz = 1:length(Zcutoffs)
    z_cutoff = Zcutoffs(zz);
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

    % titrate the equvalent threshold for delta rFP - by all simple reactions 
    % load FPA_rxns_protein_weightedFC_oriMERGE_n_100.mat
    % relFP = [relFP_f;relFP_r];
    % rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
    %               cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
    % % remove nan
    % rmInd = all(isnan(relFP),2);
    % relFP(rmInd,:) = [];
    % rowlabels(rmInd) = [];
    % relFP_wtd_local = relFP;
    % 
    % [A B] = ismember(simpleRxns,regexprep(rowlabels,'_.$','')); % if both _f and _r has rFP, _f (the previous one, is picked).
    % relFP_simpleRxns = relFP_wtd_local(B(A),:);
    % cutoffs = 0:0.01:1;
    % precisionMat = zeros(length(cutoffs),length(conditions));
    % recallMat = zeros(length(cutoffs),length(conditions));
    % OverallPrecision = zeros(length(cutoffs),1);
    % OverallRecall = zeros(length(cutoffs),1);
    % for i = 1: length(cutoffs)
    %     predictions = (relFP_simpleRxns - median(relFP_simpleRxns,2)) > cutoffs(i);
    %     TP = sum(predictions & tissueEnrichMat,1);
    %     FP = sum(predictions & ~tissueEnrichMat,1);
    %     FN = sum(~predictions & tissueEnrichMat,1);
    %     precisionMat(i,:) = TP ./ (TP + FP);
    %     recallMat(i,:) = TP ./ (TP + FN);
    %     OverallPrecision(i) = sum(TP) ./ (sum(TP) + sum(FP));
    %     OverallRecall(i) = sum(TP) ./ (sum(TP) + sum(FN));
    % end
    % precisionMat(isnan(precisionMat)) = 1; % 0/ 0;
    % figure(1)
    % hold on 
    % for i = 1:size(precisionMat,2)
    %     plot(recallMat(:,i),precisionMat(:,i),'.-');
    % end
    % legend(conditions);
    % ylabel('Precision');
    % xlabel('Recall');
    % hold off
    % 
    % figure(2)
    % hold on 
    % plot(OverallRecall, OverallPrecision,'.-');
    % ylabel('Precision');
    % xlabel('Recall');
    % plot(OverallRecall, cutoffs,'.-r');
    % hold off


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
    % figure(1)
    % hold on 
    % for i = 1:size(precisionMat,2)
    %     plot(recallMat(:,i),precisionMat(:,i),'.-');
    % end
    % legend(conditions);
    % ylabel('Precision');
    % xlabel('Recall');
    % hold off

    % figure(2)
    % hold on
    % plot(OverallRecall, OverallPrecision,'.-');
    % ylabel('Precision');
    % xlabel('Recall');
    % plot(OverallRecall, cutoffs,'.-r');
    % hold off

    % titrate the equvalent threshold for delta rFP - by unique simple reactions - raw
%     load FPA_rxns_protein_rawFC_oriMERGE_n_100.mat
%     relFP = [relFP_f;relFP_r];
%     rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
%                   cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
%     % remove nan
%     rmInd = all(isnan(relFP),2);
%     relFP(rmInd,:) = [];
%     rowlabels(rmInd) = [];
%     relFP_wtd_local = relFP;
% 
%     [A B] = ismember(simpleRxns_unique,regexprep(rowlabels,'_.$','')); % if both _f and _r has rFP, _f (the previous one, is picked).
%     relFP_simpleRxns_unique = relFP_wtd_local(B(A),:);
%     cutoffs = 0:0.01:1;
%     precisionMat_raw = zeros(length(cutoffs),length(conditions));
%     recallMat_raw = zeros(length(cutoffs),length(conditions));
%     OverallPrecision_raw = zeros(length(cutoffs),1);
%     OverallRecall_raw = zeros(length(cutoffs),1);
%     for i = 1: length(cutoffs)
%         predictions = (relFP_simpleRxns_unique - median(relFP_simpleRxns_unique,2)) >= cutoffs(i);
%         TP = sum(predictions & tissueEnrichMat_unique,1);
%         FP = sum(predictions & ~tissueEnrichMat_unique,1);
%         FN = sum(~predictions & tissueEnrichMat_unique,1);
%         precisionMat_raw(i,:) = TP ./ (TP + FP);
%         recallMat_raw(i,:) = TP ./ (TP + FN);
%         OverallPrecision_raw(i) = sum(TP) ./ (sum(TP) + sum(FP));
%         OverallRecall_raw(i) = sum(TP) ./ (sum(TP) + sum(FN));
%     end
%     precisionMat_raw(isnan(precisionMat_raw)) = 1; % 0/ 0;
    % figure(3)
    % hold on 
    % for i = 1:size(precisionMat_raw,2)
    %     plot(recallMat_raw(:,i),precisionMat_raw(:,i),'.-');
    % end
    % legend(conditions);
    % ylabel('Precision');
    % xlabel('Recall');
    % hold off

    figure(3)
    hold on
    plot(OverallRecall, OverallPrecision,'.-');
    %plot(OverallRecall_raw, OverallPrecision_raw,'.-');
    plot(OverallRecall, cutoffs,'.-');
    %plot(OverallRecall_raw, cutoffs,'.-');
    hold off


    figure(4)
    hold on
    F1 = (2 .* OverallPrecision .* OverallRecall) ./ (OverallPrecision + OverallRecall);
    %F2 = (2 .* OverallPrecision_raw .* OverallRecall_raw) ./ (OverallPrecision_raw + OverallRecall_raw);
    plot(cutoffs,F1,'.-');
    %plot(cutoffs,F2,'.-');
    hold off
end
figure(3)
ylabel('Precision');
xlabel('Recall');
legend({'PR-curve (10% TS)','cutoff-recall (10% TS)', ...
    'PR-curve (5% TS)','cutoff-recall (5% TS)',...
    'PR-curve (1% TS)','cutoff-recall (1% TS)'});
hold off

figure(4)
legend({'F-cutoff (10% TS)', 'F-cutoff (5% TS)','F-cutoff (1% TS)'});
ylabel('F-measure');
xlabel('cutoff');
hold off
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
plot(OverallRecall, OverallPrecision,'.-');
plot(OverallRecall_raw, OverallPrecision_raw,'.-');
plot(OverallRecall, cutoffs,'.-');
plot(OverallRecall_raw, cutoffs,'.-');
hold off
legend({'PR-curve (shrunk)','PR-curve (raw)','cutoff-recall (shrunk)','cutoff-recall (raw)'});
ylabel('Precison (or cutoff)')
xlabel('recall')

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
load FPA_rxns_protein_weightedFC_oriMERGE_n_100.mat
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

load FPA_rxns_protein_rawFC_oriMERGE_n_100_NAaligned.mat
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




