%% the p-value of boxplots should be changed to rank-sum test
%%
setEnvForAnalysis
addpath('PlotPub/lib')
targetExRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
targetExRxns = intersect(targetExRxns,EXinvolvedRxns);
targetExRxns = [targetExRxns;{'MI1Pt'}];
targetRxns = targetExRxns;

load('allCmp_iHumanName.mat');
allCmp_iHumanName = unique(allCmp_iHumanName);
targetRxns_DM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);

allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
targetRxns_allMetDM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);

%% load all datasets 
% protein prediction
load output/FPA_transporter_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];   
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;
load output/FPA_demand_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_DM = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_DM(rmInd) = [];
relFP_wtd_DM = relFP;
load output/pseudoDM_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
cmpName_nearest_network = cmpName_nearest;
relFP_nearest_network = relFP_nearest;
load output/FPA_allMetDemand_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM(rmInd) = [];
relFP_wtd_allMetDM = relFP;

% protein prediction - original MERGE
load output/FPA_transporter_protein_TS_common_originalFPA_originalDist_order1.5_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_oriMERGE = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];   
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_oriMERGE(rmInd) = [];
relFP_wtd_oriMERGE = relFP;
load output/FPA_demand_protein_TS_common_originalFPA_originalDist_order1.5_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_DM_oriMERGE = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_DM_oriMERGE(rmInd) = [];
relFP_wtd_DM_oriMERGE = relFP;
load output/pseudoDM_protein_TS_common_originalFPA_originalDist_order1.5_tissueNetwork.mat
cmpName_nearest_network_oriMERGE = cmpName_nearest;
relFP_nearest_network_oriMERGE = relFP_nearest;
load output/FPA_allMetDemand_protein_TS_common_originalFPA_originalDist_order1.5_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM_oriMERGE = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM_oriMERGE(rmInd) = [];
relFP_wtd_allMetDM_oriMERGE = relFP;


% RNA prediction
load output/FPA_transporter_RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_RNAcomm = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];   
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_RNAcomm(rmInd) = [];
relFP_wtd_RNAcomm = relFP;
load output/FPA_demand_RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_DM_RNAcomm = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_DM_RNAcomm(rmInd) = [];
relFP_wtd_DM_RNAcomm = relFP;
load output/pseudoDM_RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
cmpName_nearest_network_RNAcomm = cmpName_nearest;
relFP_nearest_network_RNAcomm = relFP_nearest;
load output/FPA_allMetDemand_RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM_RNAcomm = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM_RNAcomm(rmInd) = [];
relFP_wtd_allMetDM_RNAcomm = relFP;

% RNA prediction - all genes
load output/FPA_transporter_RNA_raw_all_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_RNAall = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];   
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_RNAall(rmInd) = [];
relFP_wtd_RNAall = relFP;
load output/FPA_demand_RNA_raw_all_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_DM_RNAall = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_DM_RNAall(rmInd) = [];
relFP_wtd_DM_RNAall = relFP;
load output/pseudoDM_RNA_raw_all_newFPA_weightedDist_order6_tissueNetwork.mat
cmpName_nearest_network_RNAall = cmpName_nearest;
relFP_nearest_network_RNAall = relFP_nearest;
load output/FPA_allMetDemand_RNA_raw_all_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM_RNAall = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM_RNAall(rmInd) = [];
relFP_wtd_allMetDM_RNAall = relFP;

% control - no network
load output/FPA_transporter_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_noNetwork = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_noNetwork(rmInd) = [];
relFP_wtd_noNetwork = relFP;
load output/pseudoDM_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat

%% the HMDB validation
metTissueTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','met_tissue_set_processed');
TissueAligTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','tissueAlignment');
cmpTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','metName');
allCmp = {};
validTissues = {};
for i = 1:length(TissueAligTbl.HMDBtissues)
    validTissues = [validTissues; strsplit(TissueAligTbl.HMDBtissues{i},'; ')'];
end
validTissues = unique(validTissues);
for i = 1:size(metTissueTbl,1)
    if any(strcmp(metTissueTbl.name{i},validTissues))
        cmplist = metTissueTbl{i,5};
        cmplist = strsplit(cmplist{:},'; ')';
        allCmp = union(allCmp,cmplist);
    end
end
load('metTbl.mat');
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
% match HMDB cmp to iHuman cmp
[A B] = ismember(allCmp,cmpTbl.name); 
% allCmp(~A)% failed to match, solve this manually  
allCmp_MSEAname = allCmp(A);
allCmp_HMDB = cmpTbl.hmdb_id(B(A));
[A B] = ismember(allCmp_HMDB,metTbl.HMDB);
allCmp_iHumanName = metTbl.name(B(A));
allCmp_MSEAname = allCmp_MSEAname(A);
specialMatch = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','directMatch');
allCmp_iHumanName = [allCmp_iHumanName;specialMatch.ihumanName];
allCmp_MSEAname = [allCmp_MSEAname;specialMatch.metName];


% only look at metabolites that have a transporter (the metabolite is being
% transported)
isTSP = [];
for i = 1:length(allCmp_iHumanName)
    myMet_e = {[allCmp_iHumanName{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,allCmp_iHumanName(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    candidate = intersect(myRxns_non_e,myRxns_e);
    if ~isempty(candidate)
        pass = 0;
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                pass = 1;
            end
        end
        isTSP(i) = pass;
    else
        isTSP(i) = 0;
    end

end
allCmp_iHumanName = allCmp_iHumanName(logical(isTSP));
allCmp_MSEAname = allCmp_MSEAname(logical(isTSP));

% tissue specificity matrix
measuredTissue = TissueAligTbl.name;
refMat = zeros(length(allCmp_iHumanName),length(measuredTissue));
for i = 1:length(measuredTissue)
    metList = {};
    allMets = metTissueTbl.member(ismember(metTissueTbl.name,strsplit(TissueAligTbl.HMDBtissues{i},'; ')));
    for j = 1:length(allMets)
        metList = union(metList,strsplit(allMets{j},'; '));
    end
    refMat(ismember(allCmp_MSEAname,metList),i) = 1;
end
%% delta rFP matrix - transporter FPA - protein 
relFP_sel = zeros(size(relFP_wtd,1),length(measuredTissue));
relFP_wtd_ctd = normalize(relFP_wtd,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel(:,i) = max(relFP_wtd_ctd(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

predMat = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels,'_.$',''),myRxns);

    if any(myInd)
        predMat(j,:) = max(relFP_sel(myInd,:),[],1);
    else
        predMat(j,:) = zeros(1,size(predMat,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - transporter or DM or network-FPA-nearest - protein - original MERGE
relFP_sel_oriMERGE = zeros(size(relFP_wtd_oriMERGE,1),length(measuredTissue));
relFP_wtd_ctd_oriMERGE = normalize(relFP_wtd_oriMERGE,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_oriMERGE(:,i) = max(relFP_wtd_ctd_oriMERGE(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

predMat_oriMERGE = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_oriMERGE,'_.$',''),myRxns);
    
    if any(myInd)
        predMat_oriMERGE(j,:) = max(relFP_sel_oriMERGE(myInd,:),[],1);
    else
        predMat_oriMERGE(j,:) = zeros(1,size(predMat_oriMERGE,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end
%% delta rFP matrix - transporter or DM or network-FPA-nearest - RNA common
relFP_sel_RNAcomm = zeros(size(relFP_wtd_RNAcomm,1),length(measuredTissue));
relFP_wtd_ctd_RNAcomm = normalize(relFP_wtd_RNAcomm,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_RNAcomm(:,i) = max(relFP_wtd_ctd_RNAcomm(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

predMat_RNAcomm = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_RNAcomm,'_.$',''),myRxns);
    
    if any(myInd)
        predMat_RNAcomm(j,:) = max(relFP_sel_RNAcomm(myInd,:),[],1);
    else
        predMat_RNAcomm(j,:) = zeros(1,size(predMat_RNAcomm,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - transporter or DM or network-FPA-nearest - RNA all
relFP_sel_RNAall = zeros(size(relFP_wtd_RNAall,1),length(measuredTissue));
relFP_wtd_ctd_RNAall = normalize(relFP_wtd_RNAall,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_RNAall(:,i) = max(relFP_wtd_ctd_RNAall(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

predMat_RNAall = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_RNAall,'_.$',''),myRxns);
    
    
    if any(myInd)
        predMat_RNAall(j,:) = max(relFP_sel_RNAall(myInd,:),[],1);
    else
        predMat_RNAall(j,:) = zeros(1,size(predMat_RNAall,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - transporter - no network integration
relFP_sel_noNetwork_tsp = zeros(size(relFP_wtd_noNetwork,1),length(measuredTissue));
relFP_wtd_ctd_noNetwork_tsp = normalize(relFP_wtd_noNetwork,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_noNetwork_tsp(:,i) = max(relFP_wtd_ctd_noNetwork_tsp(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

predMat_noNetwork_tsp = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_noNetwork,'_.$',''),myRxns);

    if any(myInd)
        predMat_noNetwork_tsp(j,:) = max(relFP_sel_noNetwork_tsp(myInd,:),[],1);
    else
        predMat_noNetwork_tsp(j,:) = zeros(1,size(predMat_noNetwork_tsp,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end
%% hypergeometric enrichment 
cutoffs = -1:0.01:1;
p_mat =[];
Overall_p =[];
Ncall = [];

p_mat2 =[];
Overall_p2 =[];
Ncall2 = [];

p_mat_oriMERGE =[];
Overall_p_oriMERGE =[];
Ncall_oriMERGE = [];

p_mat_RNAcomm =[];
Overall_p_RNAcomm =[];
Ncall_RNAcomm = [];

p_mat_RNAall =[];
Overall_p_RNAall =[];
Ncall_RNAall = [];

p_mat_noNetwork_tsp =[];
Overall_p_noNetwork_tsp =[];
Ncall_noNetwork_tsp = [];
for i = 1: length(cutoffs)
    predictions = predMat >= cutoffs(i) & predMat < 1.01;
    % test the narrow hypothesis 
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall(i) = sum(N);
    
    predictions = predMat_oriMERGE >= cutoffs(i)  & predMat_oriMERGE < 1.01;
    % test the narrow hypothesis 
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_oriMERGE(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_oriMERGE(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_oriMERGE(i) = sum(N);
    
    predictions = predMat_RNAcomm >= cutoffs(i) & predMat_RNAcomm < 1.01;
    % test the narrow hypothesis 
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_RNAcomm(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_RNAcomm(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_RNAcomm(i) = sum(N);
    
    predictions = predMat_RNAall >= cutoffs(i) & predMat_RNAall < 1.01;
    % test the narrow hypothesis 
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_RNAall(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_RNAall(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_RNAall(i) = sum(N);
    
    predictions = predMat_noNetwork_tsp >= cutoffs(i) & predMat_noNetwork_tsp < 1.01;
    % test the narrow hypothesis 
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_noNetwork_tsp(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_noNetwork_tsp(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_noNetwork_tsp(i) = sum(N);
    
end

figure(1)
hold on
plot(cutoffs, -log10(Overall_p),'-');
plot(cutoffs, -log10(Overall_p_noNetwork_tsp),'-');
hold off
legend({'FPA - transporter','expression - transporter '});
xlabel('cutoff')
ylabel('-log10(P value)')
xlim([-1 1]);
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [5.5, 4.5];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export('figures/benchmarkIntegrationBenefit_transporterObj.pdf');

figure(2)
hold on
plot(cutoffs, -log10(Overall_p),'-');
plot(cutoffs, -log10(Overall_p_RNAcomm),'-');
plot(cutoffs, -log10(Overall_p_RNAall),'-');
hold off
legend({'Protein FPA','RNA FPA','RNA FPA (all genes)'});
xlabel('cutoff')
ylabel('-log10(P value)')
xlim([-1 1]);
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [7, 4.5];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export('figures/benchmarkProteinVsRNA_transporterObj.svg');


figure(3)
hold on
plot(cutoffs, -log10(Overall_p_oriMERGE),'-');
plot(cutoffs, -log10(Overall_p),'-');
hold off
legend({'original FPA','improved FPA'});
xlabel('cutoff')
ylabel('-log10(P value)')
xlim([-1 1]);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [7, 4.5];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export('figures/benchmarkOriginalVsImprovedFPA_transporterObj.svg');

%% enrichment of high rFP in the tissue-specific set 
%% background delta rFP matrix - transporter FPA - protein
relFP_sel_allMetDM = zeros(size(relFP_wtd_allMetDM,1),length(measuredTissue));
relFP_wtd_ctd_allMetDM = normalize(relFP_wtd_allMetDM,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_allMetDM(:,i) = max(relFP_wtd_ctd_allMetDM(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
predMat_all = [];
qryMets = unique(metNames);
rowlabels_ID = regexprep(rowlabels,'_.$','');
S_logical = full(logical(model.S~=0));
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);

    if any(myInd)
        predMat_all = [predMat_all;max(relFP_sel(myInd,:),[],1)];
    end
end
%% background delta rFP matrix - transporter FPA - protein - original MERGE
relFP_sel_allMetDM_oriMERGE = zeros(size(relFP_wtd_allMetDM_oriMERGE,1),length(measuredTissue));
relFP_wtd_ctd_allMetDM_oriMERGE = normalize(relFP_wtd_allMetDM_oriMERGE,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_allMetDM_oriMERGE(:,i) = max(relFP_wtd_ctd_allMetDM_oriMERGE(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
predMat_all_oriMERGE = [];
rowlabels_ID_oriMERGE = regexprep(rowlabels_oriMERGE,'_.$','');
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID_oriMERGE,myRxns);
    
    if any(myInd)
        predMat_all_oriMERGE = [predMat_all_oriMERGE;max(relFP_sel_oriMERGE(myInd,:),[],1)];
    end
end
%% delta rFP matrix - transporter FPA - RNA common
relFP_sel_allMetDM_RNAcomm = zeros(size(relFP_wtd_allMetDM_RNAcomm,1),length(measuredTissue));
relFP_wtd_ctd_allMetDM_RNAcomm = normalize(relFP_wtd_allMetDM_RNAcomm,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_allMetDM_RNAcomm(:,i) = max(relFP_wtd_ctd_allMetDM_RNAcomm(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
predMat_all_RNAcomm = [];
rowlabels_ID_RNAcomm = regexprep(rowlabels_RNAcomm,'_.$','');
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID_RNAcomm,myRxns);
    
    if any(myInd)
        predMat_all_RNAcomm = [predMat_all_RNAcomm;max(relFP_sel_RNAcomm(myInd,:),[],1)];
    end
end
%% delta rFP matrix - transporter FPA - RNA all
relFP_sel_allMetDM_RNAall = zeros(size(relFP_wtd_allMetDM_RNAall,1),length(measuredTissue));
relFP_wtd_ctd_allMetDM_RNAall = normalize(relFP_wtd_allMetDM_RNAall,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_allMetDM_RNAall(:,i) = max(relFP_wtd_ctd_allMetDM_RNAall(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
predMat_all_RNAall = [];
rowlabels_ID_RNAall = regexprep(rowlabels_RNAall,'_.$','');
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID_RNAall,myRxns);
    
    if any(myInd)
        predMat_all_RNAall = [predMat_all_RNAall;max(relFP_sel_RNAall(myInd,:),[],1)];
    end
end 

%% delta rFP matrix - transporter - no network
predMat_all_noNetwork_tsp = [];
rowlabels_ID = regexprep(rowlabels_noNetwork,'_.$','');
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);
    
    DMInd = ismember(cmpName_nearest,qryMets{j});
    if any(myInd) 
        predMat_all_noNetwork_tsp = [predMat_all_noNetwork_tsp;max([relFP_sel_noNetwork_tsp(myInd,:)],[],1)];
    end
end
%% plot the box plots - benchmarkMetabolomicsPredictor_boxplot
p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP = predMat(logical(refMat));
boxData = [boxData;enrichedMet_rFP;predMat_all(:)];
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1)*size(predMat_all,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{met+rxn}';'all metabolites^{met+rxn}'}];
p = ranksum(enrichedMet_rFP,predMat_all(:),  'Tail','right');

%% plot the box plots - benchmarkIntegration_boxplot
p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP_noNetwork_tsp = predMat_noNetwork_tsp(logical(refMat));
boxData = [boxData;enrichedMet_rFP_noNetwork_tsp;predMat_all_noNetwork_tsp(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_noNetwork_tsp'},length(enrichedMet_rFP_noNetwork_tsp),1);repmat({'all_noNetwork_tsp'},size(predMat_all_noNetwork_tsp,1)*size(predMat_all_noNetwork_tsp,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{tranporter expression}';'all metabolites^{tranporter expression}'}];
p(1) = ranksum(enrichedMet_rFP_noNetwork_tsp,predMat_all_noNetwork_tsp(:),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat));
boxData = [boxData;enrichedMet_rFP;predMat_all(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1)*size(predMat_all,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{FPA}';'all metabolites^{FPA}'}];
p(2) = ranksum(enrichedMet_rFP,predMat_all(:),  'Tail','right');

figure('units','inch','position',[0,0,9,8])
boxplot(boxData, boxLabel,'Labels',xLabels,'Colors','bkbk','Symbol','ro');
ylim([-1.1 1.1])
ylabel('\DeltarFP');
xlabel('');
set(gca, 'TickLabelInterpreter', 'tex');
text(1.3,0.75,['p < ',num2str(p(1),2)],'fontsize', 12,'BackgroundColor','w')
text(3.3,0.75,['p < ',num2str(p(2),2)],'fontsize', 12,'BackgroundColor','w')
set(gca, 'fontsize', 14);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xtickangle(30)
S = hgexport('readstyle','default_sci');
style.Format = 'pdf';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
saveas(gca,'figures/benchmarkIntegrationBenefit_boxplot_transporterObj.pdf');

%% plot the box plots - benchmarkProteinVsRNA_boxplot
p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP = predMat(logical(refMat));
boxData = [boxData;enrichedMet_rFP;predMat_all(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1)*size(predMat_all,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{protein}';'all metabolites^{protein}'}];
p(1) = ranksum(enrichedMet_rFP,predMat_all(:),  'Tail','right');

enrichedMet_rFP_RNAcomm = predMat_RNAcomm(logical(refMat));
boxData = [boxData;enrichedMet_rFP_RNAcomm;predMat_all_RNAcomm(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_RNAcomm'},length(enrichedMet_rFP_RNAcomm),1);repmat({'all_RNAcomm'},size(predMat_all_RNAcomm,1)*size(predMat_all_RNAcomm,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{RNA}';'all metabolites^{RNA}'}];
p(2) = ranksum(enrichedMet_rFP_RNAcomm,predMat_all_RNAcomm(:),  'Tail','right');

enrichedMet_rFP_RNAall = predMat_RNAall(logical(refMat));
boxData = [boxData;enrichedMet_rFP_RNAall;predMat_all_RNAall(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_RNAall'},length(enrichedMet_rFP_RNAall),1);repmat({'all_RNAall'},size(predMat_all_RNAall,1)*size(predMat_all_RNAall,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{RNA (all genes)}';'all metabolites^{RNA (all genes)}'}];
p(3) = ranksum(enrichedMet_rFP_RNAall,predMat_all_RNAall(:),  'Tail','right');


figure('units','inch','position',[0,0,10,10])
boxplot(boxData, boxLabel,'Labels',xLabels,'Colors','bkbk','Symbol','ro');
ylim([-1.1 1.1])
ylabel('\DeltarFP');
xlabel('');
set(gca, 'TickLabelInterpreter', 'tex');
text(1.3,0.75,['p < ',num2str(p(1),2)],'fontsize', 12,'BackgroundColor','w')
text(3.3,0.75,['p < ',num2str(p(2),2)],'fontsize', 12,'BackgroundColor','w')
text(5.3,0.75,['p < ',num2str(p(3),2)],'fontsize', 12,'BackgroundColor','w')
set(gca, 'fontsize', 14);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xtickangle(30)
S = hgexport('readstyle','default_sci');
style.Format = 'svg';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
saveas(gca,'figures/benchmarkProteinVsRNA_boxplot_transporterObj.pdf');

%% plot the box plots - benchmarkOriginalVsNewFPA_boxplot
p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP_oriMERGE = predMat_oriMERGE(logical(refMat));
boxData = [boxData;enrichedMet_rFP_oriMERGE;predMat_all_oriMERGE(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched_oriMERGE'},length(enrichedMet_rFP_oriMERGE),1);repmat({'all_oriMERGE'},size(predMat_all_oriMERGE,1)*size(predMat_all_oriMERGE,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{original FPA}';'all metabolites^{original FPA}'}];
p(1) = ranksum(enrichedMet_rFP_oriMERGE,predMat_all_oriMERGE(:),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat));
boxData = [boxData;enrichedMet_rFP;predMat_all(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1)*size(predMat_all,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{new FPA}';'all metabolites^{new FPA}'}];
p(2) = ranksum(enrichedMet_rFP,predMat_all(:),  'Tail','right');

figure('units','inch','position',[0,0,8,8])
boxplot(boxData, boxLabel,'Labels',xLabels,'Colors','bkbk','Symbol','ro');
ylim([-1.1 1.1])
ylabel('\DeltarFP');
xlabel('');
set(gca, 'TickLabelInterpreter', 'tex');
text(1.3,0.75,['p < ',num2str(p(1),2)],'fontsize', 12)
text(3.3,0.75,['p < ',num2str(p(2),2)],'fontsize', 12)
set(gca, 'fontsize', 14);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xtickangle(30)
S = hgexport('readstyle','default_sci');
style.Format = 'svg';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
saveas(gca,'figures/benchmarkOriginalVsNewFPA_boxplot_transporterObj.pdf');
%% box plot of individual tissues 
% i = 1, Adrenal Gland' is predicted by expression of related rxns
% i = 5 skin shows the benefit of integration, other examples are: lung (i
% = 9), prostate (i = 13), testis i = 14
% in many tissues like intestine and liver, the FPA IS NOT outperforming
% tsp + rxn expression
i = 17;

measuredTissue(i)
sum(logical(refMat(:,i)))

p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP_noNetwork_tsp = predMat_noNetwork_tsp(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP_noNetwork_tsp;predMat_all_noNetwork_tsp(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_noNetwork_tsp'},length(enrichedMet_rFP_noNetwork_tsp),1);repmat({'all_noNetwork_tsp'},size(predMat_all_noNetwork_tsp,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{tranporter expression}';'all metabolites^{tranporter expression}'}];
p(1) = ranksum(enrichedMet_rFP_noNetwork_tsp,predMat_all_noNetwork_tsp(:,i),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP;predMat_all(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{FPA}';'all metabolites^{FPA}'}];
p(2) = ranksum(enrichedMet_rFP,predMat_all(:,i),  'Tail','right');


figure('units','inch','position',[0,0,8,8])
boxplot(boxData, boxLabel,'Labels',xLabels,'Colors','bkbk','Symbol','ro');
ylim([-1.1 1.1])
ylabel('\DeltarFP');
xlabel('');
set(gca, 'TickLabelInterpreter', 'tex');
text(1.3,0.75,['p < ',num2str(p(1),2)],'fontsize', 12,'BackgroundColor','w')
text(3.3,0.75,['p < ',num2str(p(2),2)],'fontsize', 12,'BackgroundColor','w')
set(gca, 'fontsize', 14);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xtickangle(30)
S = hgexport('readstyle','default_sci');
style.Format = 'svg';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
%% save some desired ones
saveas(gca,['figures/benchmarkIntegrationBenefit_boxplot_individualTissue_',measuredTissue{i},'_transporterObj.pdf']);

%% box plot of individual tissues 
% i = 1, Adrenal Gland' RNA outperforms protein; also i = 2, Artery (but
% only 2 mets)
% i = 5, skin ; protein is better; i = 8 intestine 
i = 8;

measuredTissue(i)
sum(logical(refMat(:,i)))

p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP = predMat(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP;predMat_all(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{protein}';'all metabolites^{protein}'}];
p(1) = ranksum(enrichedMet_rFP,predMat_all(:,i),  'Tail','right');

enrichedMet_rFP_RNAcomm = predMat_RNAcomm(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP_RNAcomm;predMat_all_RNAcomm(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_RNAcomm'},length(enrichedMet_rFP_RNAcomm),1);repmat({'all_RNAcomm'},size(predMat_all_RNAcomm,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{RNA}';'all metabolites^{RNA}'}];
p(2) = ranksum(enrichedMet_rFP_RNAcomm,predMat_all_RNAcomm(:,i),  'Tail','right');

enrichedMet_rFP_RNAall = predMat_RNAall(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP_RNAall;predMat_all_RNAall(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_RNAall'},length(enrichedMet_rFP_RNAall),1);repmat({'all_RNAall'},size(predMat_all_RNAall,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{RNA (all genes)}';'all metabolites^{RNA (all genes)}'}];
p(3) = ranksum(enrichedMet_rFP_RNAall,predMat_all_RNAall(:,i),  'Tail','right');


figure('units','inch','position',[0,0,8,8])
boxplot(boxData, boxLabel,'Labels',xLabels,'Colors','bkbk','Symbol','ro');
ylim([-1.1 1.1])
ylabel('\DeltarFP');
xlabel('');
set(gca, 'TickLabelInterpreter', 'tex');
text(1.3,0.75,['p < ',num2str(p(1),2)],'fontsize', 12,'BackgroundColor','w')
text(3.3,0.75,['p < ',num2str(p(2),2)],'fontsize', 12,'BackgroundColor','w')
text(5.3,0.75,['p < ',num2str(p(3),2)],'fontsize', 12,'BackgroundColor','w')
set(gca, 'fontsize', 14);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xtickangle(30)
S = hgexport('readstyle','default_sci');
style.Format = 'svg';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
%% save some desired ones
saveas(gca,['figures/benchmarkProteinVsRNA_boxplot_individualTissue_',measuredTissue{i},'_transporterObj.svg']);