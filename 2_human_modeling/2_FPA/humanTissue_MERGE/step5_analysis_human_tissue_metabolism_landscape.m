%% About
% generate the predictions of tissue-enriched metabolic fluxes and save the
% tables for further visualization
%% load environmnet
setEnvForAnalysis
% dataset
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

%% filter the predictions by the delta FC cutoff of 0.2
relFP_filtered = relFP_wtd - median(relFP_wtd,2);
rowlabels_filtered = rowlabels;
relFP_filtered_RNA = relFP_RNA - median(relFP_RNA,2);
keep = any(abs(relFP_filtered) > 0.2 | abs(relFP_filtered_RNA) > 0.2,2);
relFP_filtered = relFP_filtered(keep,:);
rowlabels_filtered = rowlabels_filtered(keep);
relFP_filtered_RNA = relFP_filtered_RNA(keep,:);

length(rowlabels_filtered)
%% filter the predictions by the consistency with iMAT++ OFD (see Yilmaz et al 2020 for more info about the MERGE pipelien)
% the main purpose of this consistency check is to narrow down the
% prediction to a specific direction of the reaction, when the reaction is
% reversibe
% load the OFD
lib = 'protein';
for i = 1:length(conditions)
    load(['./../../1_iMAT++/output/humanModel/',lib,'/OFD/',conditions{i},'.mat']);
    eval([conditions{i},'=myCSM;']);
end
% check consistency with OFD
IsInOFD = false(length(rowlabels_filtered),1);
for j = 1:length(rowlabels_filtered)
    % we consider the reaction is consistent as long as one of the
    % enriched-tissue (delta rFP>0.2) is in its OFD
    highestSites = conditions((relFP_filtered(j,:) > 0.2));
    for i = 1:length(highestSites)
        eval(['myCSM=',highestSites{i},';']);
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
        else % we now only look at internal reactions, return error
            error('wrong input!');
        end
    end
end
IsInOFD_pro = IsInOFD;

% the same consistency check for RNA
lib = 'RNA';
for i = 1:length(conditions)
    load(['./../../1_iMAT++/output/humanModel/',lib,'/OFD/',conditions{i},'.mat']);
    eval([conditions{i},'=myCSM;']);
end
IsInOFD = false(length(rowlabels_filtered),1);
for j = 1:length(rowlabels_filtered)
    highestSites = conditions((relFP_filtered_RNA(j,:) > 0.2));
    for i = 1:length(highestSites)
        eval(['myCSM=',highestSites{i},';']);
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
        else % we now only look at internal reactions, return error
            error('wrong input!');
        end
    end
end
IsInOFD_RNA = IsInOFD;

% sum(IsInOFD_pro)
% sum(IsInOFD_RNA)
% we filter out the predictions that are neither consistent with
% protein-based nor RNA-based OFD
relFP_filtered = relFP_filtered(IsInOFD_pro | IsInOFD_RNA,:);
relFP_filtered_RNA = relFP_filtered_RNA(IsInOFD_pro | IsInOFD_RNA,:);
rowlabels_filtered = rowlabels_filtered(IsInOFD_pro | IsInOFD_RNA,:);

%% plot relative FP clustergram
distMethod_col = 'cosine';
distMethod_row = 'cosine';

IDs = regexprep(conditions,'_','-');
IDs = [cellfun(@(x) [x,' - P'],IDs,'UniformOutput',0),cellfun(@(x) [x,' - R'],IDs,'UniformOutput',0)];

relFP_filtered2 = [relFP_filtered,relFP_filtered_RNA];

cgo=clustergram(relFP_filtered2,'RowLabels',rowlabels_filtered,'ColumnLabels',IDs,'RowPDist',distMethod_row,'ColumnPDist',distMethod_col);
% set the plot attributes
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
set(h, 'FontSize', 12)
% save figure
S = hgexport('readstyle','default_sci');
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
% save
set(h, 'FontSize', 7)
exportgraphics(gcf,'figures/clustergram_rxn_FPA.tiff','ContentType','vector')
exportgraphics(gcf,'figures/clustergram_rxn_FPA.pdf','ContentType','vector')
%% inspect the clusters and label information (this is interactive)
ROI = mygroup.RowNodeNames;
ROI = regexprep(ROI,'_.$','');
[ROI_annotated, ROI_enrichment] = annotateRxnSet(ROI,model);

%% visualize subsystems of enriched rxns for each tissue 
subsys = unique(model.subSystems);
subsys = setdiff(subsys,{''});
countMat = zeros(length(subsys),length(conditions));
n_pop = length(targetRxns);
popSubSys = model.subSystems(ismember(model.rxns,targetRxns));
for i = 1:length(conditions)
    tissue = conditions{i};
    ROI = rowlabels_filtered(relFP_filtered2(:,find(strcmp(conditions,tissue))) >= 0.2); % we consider a reaction tissue-enriched as long as it is enriched in one dataset (protein or rna) 
    ROI = regexprep(ROI,'_.$','');
    % calculate percentage
    subsystemInfo = {};
    subsystemInfo(:,1) = unique(model.subSystems(ismember(model.rxns,ROI)));
    rxnInd = ismember(model.rxns,ROI);
    count = [];
    for j = 1:length(subsystemInfo(:,1))
        mysys = subsystemInfo{j,1};
        n_obs = sum(strcmp(model.subSystems(rxnInd),mysys));
        n_trueInPop = sum(ismember(popSubSys ,mysys));
        count(j) = n_obs;
    end
    [A B] = ismember(subsys,subsystemInfo(:,1));
    countMat(A,i) = count(B(A));
end
% remove housekeeping (no enrichment in any tissue)
housekeepPathways = subsys(max(countMat,[],2) < 1)
subsys_enri = subsys(max(countMat,[],2) > 0);% remove housekeeping
countMat_enri = countMat(max(countMat,[],2) > 0,:);% remove housekeeping
countMat_enri = countMat_enri ./ max(countMat_enri,[],2); % row-wise normalize

% make the publishable figure in R
t = array2table(countMat_enri);
t.Properties.RowNames = subsys_enri;
writetable(t,'output/subsys_annotation_for_enriched_rxns_heatmap.csv','WriteRowNames',1);

boundaries =  cell2table(conditions');
writetable(boundaries,'output/subsys_annotation_for_enriched_rxns_colnames.csv');



%% determination of the delta rFP cutoff - why we choose to use 0.2 as cutoff
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

    %%titrate the equvalent threshold for delta rFP - by unique simple reactions 
    load output/FPA_rxn_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
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
    precisionMat(isnan(precisionMat)) = 1; % no true positive, we skip
    
    figure(3)
    hold on
    plot(OverallRecall, OverallPrecision,'.-');
    plot(OverallRecall, cutoffs,'.-');
    hold off


    figure(4)
    hold on
    F1 = (2 .* OverallPrecision .* OverallRecall) ./ (OverallPrecision + OverallRecall);% F1 metric of performance
    plot(cutoffs,F1,'.-');
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

