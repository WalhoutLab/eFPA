%% About 
% systematically benchmark FPA (the principle of local expression
% dictating metabolic flux/functions) and compare protein data with RNA
% data as the input, based on HMDB reference tissue-enriched metabolite
% dataset
%% set up environments
setEnvForAnalysis
addpath('./PlotPub/lib')
targetExRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
targetExRxns = intersect(targetExRxns,EXinvolvedRxns);
targetExRxns = [targetExRxns;{'MI1Pt'}];
targetRxns = targetExRxns; % in this script, the targetrxns by default refers to the transporter objectives
% load all HMDB reference cmps
load('allCmp_iHumanName.mat');
allCmp_iHumanName = unique(allCmp_iHumanName);
targetRxns_DM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);
% all metabolites objectives
allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
targetRxns_allMetDM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);
%% load FPA modeling result datasets 
% protein prediction
% transporters
load output/FPA_transporter_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];   
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;
% demand and sinks for HMDB reference 
load output/FPA_demand_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_DM = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_DM(rmInd) = [];
relFP_wtd_DM = relFP;
% associated internal reactions
load output/pseudoDM_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
cmpName_nearest_network = cmpName_nearest;
relFP_nearest_network = relFP_nearest;
% demand and sink for all mets
load output/FPA_allMetDemand_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM(rmInd) = [];
relFP_wtd_allMetDM = relFP;

% RNA prediction (by common genes only)
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

% RNA prediction (by all genes)
load output/FPA_transporter_RNA_TS_all_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_RNAall = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];   
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_RNAall(rmInd) = [];
relFP_wtd_RNAall = relFP;
load output/FPA_demand_RNA_TS_all_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_DM_RNAall = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_DM_RNAall(rmInd) = [];
relFP_wtd_DM_RNAall = relFP;
load output/pseudoDM_RNA_TS_all_newFPA_weightedDist_order6_tissueNetwork.mat
cmpName_nearest_network_RNAall = cmpName_nearest;
relFP_nearest_network_RNAall = relFP_nearest;
load output/FPA_allMetDemand_RNA_TS_all_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM_RNAall = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM_RNAall(rmInd) = [];
relFP_wtd_allMetDM_RNAall = relFP;

% control - target expression only, no network
load output/FPA_transporter_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_noNetwork = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_noNetwork(rmInd) = [];
relFP_wtd_noNetwork = relFP;
load output/pseudoDM_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
cmpName_nearest_network_noNetwork = cmpName_nearest;
relFP_nearest_network_noNetwork = relFP_nearest;
%% load the HMDB reference tissue-enriched metabolite set 
% load HMDB data
metTissueTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','met_tissue_set_processed');
TissueAligTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','tissueAlignment');
cmpTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','metName');
% extract the overlapped tissues and compound names
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
allCmp_MSEAname = allCmp(A);
allCmp_HMDB = cmpTbl.hmdb_id(B(A));
[A B] = ismember(allCmp_HMDB,metTbl.HMDB);
allCmp_iHumanName = metTbl.name(B(A));
allCmp_MSEAname = allCmp_MSEAname(A);
specialMatch = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','directMatch');
allCmp_iHumanName = [allCmp_iHumanName;specialMatch.ihumanName];
allCmp_MSEAname = [allCmp_MSEAname;specialMatch.metName];

% finally, make the tissue-enriched metabolite matrix
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
%% make the delta rFP matrix - protein data - general met predictor (tsp, dm/snk, asso. internal.) 
% transporter rFP
relFP_sel = zeros(size(relFP_wtd,1),length(measuredTissue));
relFP_wtd_ctd = normalize(relFP_wtd,2,'center','median');
% we first merge the similar tissues by taking the maximal delta rFP
for i = 1:length(measuredTissue)
    relFP_sel(:,i) = max(relFP_wtd_ctd(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
% same for demand rFP
relFP_sel_DM = zeros(size(relFP_wtd_DM,1),length(measuredTissue));
relFP_wtd_ctd_DM = normalize(relFP_wtd_DM,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_DM(:,i) = max(relFP_wtd_ctd_DM(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
% and same for associated internal reactions
relFP_sel_nearest = zeros(size(relFP_nearest_network,1),length(measuredTissue));
relFP_wtd_ctd_nearest = relFP_nearest_network;
for i = 1:length(measuredTissue)
    relFP_sel_nearest(:,i) = max(relFP_wtd_ctd_nearest(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

% next, we take the maximum of differnet objectives
predMat = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels,'_.$',''),myRxns);
    DMrxn = {['NewMet_',allCmp_iHumanName{j},'_f'],['NewMet_',allCmp_iHumanName{j},'_r']};
    DMInd = ismember(rowlabels_DM,DMrxn);
    NearestInd = ismember(cmpName_nearest_network,allCmp_iHumanName{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat(j,:) = max([relFP_sel(myInd,:);relFP_sel_DM(DMInd,:);relFP_sel_nearest(NearestInd,:)],[],1);
    else
        predMat(j,:) = zeros(1,size(predMat,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - protein data - met-centric predictor (transporter or DM)
% for benchmakring the predictor, we tested the prediction when excluding 
% information from internal associated rxns
predMat2 = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels,'_.$',''),myRxns);
    DMrxn = {['NewMet_',allCmp_iHumanName{j},'_f'],['NewMet_',allCmp_iHumanName{j},'_r']};
    DMInd = ismember(rowlabels_DM,DMrxn);
        
    if any(myInd) || any(DMInd)
        predMat2(j,:) = max([relFP_sel(myInd,:);relFP_sel_DM(DMInd,:)],[],1);
    else
        predMat2(j,:) = zeros(1,size(predMat2,2));
        fprintf('%s is not predicted without internal associated rxns\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - RNA common - general predictor (transporter or DM or network-FPA-nearest)
% merge similar tissues
relFP_sel_RNAcomm = zeros(size(relFP_wtd_RNAcomm,1),length(measuredTissue));
relFP_wtd_ctd_RNAcomm = normalize(relFP_wtd_RNAcomm,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_RNAcomm(:,i) = max(relFP_wtd_ctd_RNAcomm(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_DM_RNAcomm = zeros(size(relFP_wtd_DM_RNAcomm,1),length(measuredTissue));
relFP_wtd_ctd_DM_RNAcomm = normalize(relFP_wtd_DM_RNAcomm,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_DM_RNAcomm(:,i) = max(relFP_wtd_ctd_DM_RNAcomm(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_nearest_RNAcomm = zeros(size(relFP_nearest_network_RNAcomm,1),length(measuredTissue));
relFP_wtd_ctd_nearest_RNAcomm = relFP_nearest_network_RNAcomm;
for i = 1:length(measuredTissue)
    relFP_sel_nearest_RNAcomm(:,i) = max(relFP_wtd_ctd_nearest_RNAcomm(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
% take maximum among the three objectives
predMat_RNAcomm = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_RNAcomm,'_.$',''),myRxns);
    DMrxn = {['NewMet_',allCmp_iHumanName{j},'_f'],['NewMet_',allCmp_iHumanName{j},'_r']};
    DMInd = ismember(rowlabels_DM_RNAcomm,DMrxn);
    NearestInd = ismember(cmpName_nearest_network_RNAcomm,allCmp_iHumanName{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_RNAcomm(j,:) = max([relFP_sel_RNAcomm(myInd,:);relFP_sel_DM_RNAcomm(DMInd,:);relFP_sel_nearest_RNAcomm(NearestInd,:)],[],1);
    else
        predMat_RNAcomm(j,:) = zeros(1,size(predMat_RNAcomm,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - RNA all genes - general predictor (transporter or DM or network-FPA-nearest)
% merge similar tissues
relFP_sel_RNAall = zeros(size(relFP_wtd_RNAall,1),length(measuredTissue));
relFP_wtd_ctd_RNAall = normalize(relFP_wtd_RNAall,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_RNAall(:,i) = max(relFP_wtd_ctd_RNAall(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_DM_RNAall = zeros(size(relFP_wtd_DM_RNAall,1),length(measuredTissue));
relFP_wtd_ctd_DM_RNAall = normalize(relFP_wtd_DM_RNAall,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_DM_RNAall(:,i) = max(relFP_wtd_ctd_DM_RNAall(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_nearest_RNAall = zeros(size(relFP_nearest_network_RNAall,1),length(measuredTissue));
relFP_wtd_ctd_nearest_RNAall = relFP_nearest_network_RNAall;
for i = 1:length(measuredTissue)
    relFP_sel_nearest_RNAall(:,i) = max(relFP_wtd_ctd_nearest_RNAall(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
% take maximum of the three objectives
predMat_RNAall = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_RNAall,'_.$',''),myRxns);
    DMrxn = {['NewMet_',allCmp_iHumanName{j},'_f'],['NewMet_',allCmp_iHumanName{j},'_r']};
    DMInd = ismember(rowlabels_DM_RNAall,DMrxn);
    NearestInd = ismember(cmpName_nearest_network_RNAall,allCmp_iHumanName{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_RNAall(j,:) = max([relFP_sel_RNAall(myInd,:);relFP_sel_DM_RNAall(DMInd,:);relFP_sel_nearest_RNAall(NearestInd,:)],[],1);
    else
        predMat_RNAall(j,:) = zeros(1,size(predMat_RNAall,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - no network integration - general predictor (transporter or asso. int. rxns)
% first merge similar tissues
relFP_sel_noNetwork = zeros(size(relFP_wtd_noNetwork,1),length(measuredTissue));
relFP_wtd_ctd_noNetwork = normalize(relFP_wtd_noNetwork,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_noNetwork(:,i) = max(relFP_wtd_ctd_noNetwork(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_DM_noNetwork = zeros(size(relFP_nearest_noNetwork,1),length(measuredTissue));
relFP_wtd_ctd_DM_NoNetwork = relFP_nearest_noNetwork;% DM/SNK is artificial reaction for FBA-based analysis, so here we just use the delta rFP based on all associated internal reactions' expression 
for i = 1:length(measuredTissue)
    relFP_sel_DM_noNetwork(:,i) = max(relFP_wtd_ctd_DM_NoNetwork(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
% take maximum 
predMat_noNetwork = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_noNetwork,'_.$',''),myRxns);
    DMInd = ismember(cmpName_nearest_noNetwork,allCmp_iHumanName{j});
    if any(myInd) || any(DMInd)
        predMat_noNetwork(j,:) = max([relFP_sel_noNetwork(myInd,:);relFP_sel_DM_noNetwork(DMInd,:)],[],1);
    else
        predMat_noNetwork(j,:) = zeros(1,size(predMat_noNetwork,2));
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end
%% delta rFP matrix - no network integration - transporter expression only
% merge the tissues
relFP_sel_noNetwork_tsp = zeros(size(relFP_wtd_noNetwork,1),length(measuredTissue));
relFP_wtd_ctd_noNetwork_tsp = normalize(relFP_wtd_noNetwork,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_noNetwork_tsp(:,i) = max(relFP_wtd_ctd_noNetwork_tsp(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
% take maximum among all associated transporters
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
%% assessing the prediction: hypergeometric enrichment of the overlap 
cutoffs = -1:0.01:1;
% protein - common genes - this is considered as default 
p_mat =[];
Overall_p =[];
Ncall = [];
% protein - common genes - without asso. int. rxns
p_mat2 =[];
Overall_p2 =[];
Ncall2 = [];

p_mat_RNAcomm =[];
Overall_p_RNAcomm =[];
Ncall_RNAcomm = [];

p_mat_RNAall =[];
Overall_p_RNAall =[];
Ncall_RNAall = [];

p_mat_noNetwork =[];
Overall_p_noNetwork =[];
Ncall_noNetwork = [];

p_mat_noNetwork_tsp =[];
Overall_p_noNetwork_tsp =[];
Ncall_noNetwork_tsp = [];
for i = 1: length(cutoffs)
    predictions = predMat >= cutoffs(i) & predMat < 1.01;% if it is larger than 1, we basically consider it numerical error
    % although the numerical-error related predictions are likely to be a
    % real high delta rFP, we exclude it to be conservative. Since human1
    % model is too large, few numeric errors in a large-scale computation
    % is hard to avoid. 
    
    % test the enrichment
    x = sum(predictions & refMat,1); % overlap number 
    M = size(predictions,1); % total number of cases
    K = sum(refMat,1); % known tissue-enriched metabolite labels
    N = sum(predictions,1); % predicted enriched labels 
    p_mat(i,:) = hygecdf(x-1,M,K,N,'upper');% tissue-by-tissue enrichment
    Overall_p(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall(i) = sum(N);
    
    % apply the same calculation to other predictions
    predictions = predMat2 >= cutoffs(i) & predMat2 < 1.01;
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat2(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p2(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall2(i) = sum(N);
    
    predictions = predMat_RNAcomm >= cutoffs(i) & predMat_RNAcomm < 1.01;
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_RNAcomm(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_RNAcomm(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_RNAcomm(i) = sum(N);
    
    predictions = predMat_RNAall >= cutoffs(i) & predMat_RNAall < 1.01;
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_RNAall(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_RNAall(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_RNAall(i) = sum(N);
    
    predictions = predMat_noNetwork >= cutoffs(i) & predMat_noNetwork < 1.01;
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_noNetwork(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_noNetwork(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_noNetwork(i) = sum(N);
    
    predictions = predMat_noNetwork_tsp >= cutoffs(i) & predMat_noNetwork_tsp < 1.01;
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
plot(cutoffs, -log10(Overall_p2),'-');
plot(cutoffs, -log10(Overall_p),'-');
hold off
legend({'FPA - exchaning','FPA - exchaning+internal'});
xlabel('cutoff')
ylabel('-log10(P value)')
xlim([-1 1]);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3.5, 2.25];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.TickDir = 'out';
plt.AxisLineWidth = 1;
plt.export('figures/benchmarkMetabolomicsPredictor.pdf');

figure(2)
hold on
plot(cutoffs, -log10(Overall_p),'-');
plot(cutoffs, -log10(Overall_p_noNetwork_tsp),'-');
plot(cutoffs, -log10(Overall_p_noNetwork),'-');
hold off
legend({'FPA','Transporter expression',sprintf('Expression of transporter \nand associated rxn(s)')});
xlabel('cutoff')
ylabel('-log10(P value)')
xlim([-1 1]);
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [3.5, 2.25];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.TickDir = 'out';
plt.AxisLineWidth = 1;
plt.export('figures/benchmarkIntegrationBenefit.pdf');

figure(3)
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
plt.BoxDim = [3.5, 2.25];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.AxisLineWidth = 1;
plt.TickDir = 'out';
plt.export('figures/benchmarkProteinVsRNA.pdf');

%% assessing the delta rFP directly 
% i.e. by enrichment of high rFP in the tissue-specific set 
%% load the background delta rFP - protein comm - general predictor (transporter or DM or nearest networkFPA)
% same processing as above, just apply to all the metabolites in human 1
% model
% first merge tissues for the demand/snk objective FPA predictions of all
% metabolites
relFP_sel_allMetDM = zeros(size(relFP_wtd_allMetDM,1),length(measuredTissue));
relFP_wtd_ctd_allMetDM = normalize(relFP_wtd_allMetDM,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_allMetDM(:,i) = max(relFP_wtd_ctd_allMetDM(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
predMat_all = [];
qryMets = unique(metNames);% all unique mets in the model
rowlabels_ID = regexprep(rowlabels,'_.$','');% rowlabels is the rxn ID for the FPA of all transporters
S_logical = full(logical(model.S~=0));
% take maximum among the three types of objectives
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);
    DMrxn = {['NewMet_',qryMets{j},'_f'],['NewMet_',qryMets{j},'_r']};
    DMInd = ismember(rowlabels_allMetDM,DMrxn);
    NearestInd = ismember(cmpName_nearest_network,qryMets{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_all = [predMat_all;max([relFP_sel(myInd,:);relFP_sel_allMetDM(DMInd,:);relFP_sel_nearest(NearestInd,:)],[],1)];
    end
end

%% background delta rFP matrix - protein comm - no int. asso. rxn (transporter or DM)
predMat_all2 = [];
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);
    DMrxn = {['NewMet_',qryMets{j},'_f'],['NewMet_',qryMets{j},'_r']};
    DMInd = ismember(rowlabels_allMetDM,DMrxn);
        
    if any(myInd) || any(DMInd)
        predMat_all2 = [predMat_all2;max([relFP_sel(myInd,:);relFP_sel_allMetDM(DMInd,:)],[],1)];
    end
end
%% background delta rFP matrix - RNA common - general predictor (transporter or DM or network-FPA-nearest)
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
    DMrxn = {['NewMet_',qryMets{j},'_f'],['NewMet_',qryMets{j},'_r']};
    DMInd = ismember(rowlabels_allMetDM_RNAcomm,DMrxn);
    NearestInd = ismember(cmpName_nearest_network_RNAcomm,qryMets{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_all_RNAcomm = [predMat_all_RNAcomm;max([relFP_sel_RNAcomm(myInd,:);relFP_sel_allMetDM_RNAcomm(DMInd,:);relFP_sel_nearest_RNAcomm(NearestInd,:)],[],1)];
    end
end
%% background delta rFP matrix - RNA all - general predictor (transporter or DM or network-FPA-nearest)
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
    DMrxn = {['NewMet_',qryMets{j},'_f'],['NewMet_',qryMets{j},'_r']};
    DMInd = ismember(rowlabels_allMetDM_RNAall,DMrxn);
    NearestInd = ismember(cmpName_nearest_network_RNAall,qryMets{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_all_RNAall = [predMat_all_RNAall;max([relFP_sel_RNAall(myInd,:);relFP_sel_allMetDM_RNAall(DMInd,:);relFP_sel_nearest_RNAall(NearestInd,:)],[],1)];
    end
end
%% background delta rFP matrix - no network - general predictor (transporter or asso. int. rxns)
predMat_all_noNetwork = [];
rowlabels_ID = regexprep(rowlabels_noNetwork,'_.$','');
S_logical = full(logical(model.S~=0));
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);
    DMInd = ismember(cmpName_nearest_noNetwork,qryMets{j});
    if any(myInd) || any(DMInd)
        predMat_all_noNetwork = [predMat_all_noNetwork;max([relFP_sel_noNetwork(myInd,:);relFP_sel_DM_noNetwork(DMInd,:)],[],1)];
    end
end
%% background delta rFP matrix - no network - transporter only
predMat_all_noNetwork_tsp = [];
rowlabels_ID = regexprep(rowlabels_noNetwork,'_.$','');
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);
    if any(myInd) 
        predMat_all_noNetwork_tsp = [predMat_all_noNetwork_tsp;max([relFP_sel_noNetwork(myInd,:)],[],1)];
    end
end
%% plot the box plots
p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP2 = predMat2(logical(refMat));
boxData = [boxData;enrichedMet_rFP2;predMat_all2(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched2'},length(enrichedMet_rFP2),1);repmat({'all2'},size(predMat_all2,1)*size(predMat_all2,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{met}';'all metabolites^{met}'}];
p(1) = ranksum(enrichedMet_rFP2,predMat_all2(:),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat));
boxData = [boxData;enrichedMet_rFP;predMat_all(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1)*size(predMat_all,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{met+rxn}';'all metabolites^{met+rxn}'}];
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
saveas(gca,'figures/benchmarkMetabolomicsPredictor_boxplot.pdf');
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

enrichedMet_rFP_noNetwork = predMat_noNetwork(logical(refMat));
boxData = [boxData;enrichedMet_rFP_noNetwork;predMat_all_noNetwork(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_noNetwork'},length(enrichedMet_rFP_noNetwork),1);repmat({'all_noNetwork'},size(predMat_all_noNetwork,1)*size(predMat_all_noNetwork,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{tranporter+rxn expression}';'all metabolites^{tranporter+rxn expression}'}];
p(2) = ranksum(enrichedMet_rFP_noNetwork,predMat_all_noNetwork(:),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat));
boxData = [boxData;enrichedMet_rFP;predMat_all(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1)*size(predMat_all,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{FPA}';'all metabolites^{FPA}'}];
p(3) = ranksum(enrichedMet_rFP,predMat_all(:),  'Tail','right');

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
saveas(gca,'figures/benchmarkIntegrationBenefit_boxplot.pdf');

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
saveas(gca,'figures/benchmarkProteinVsRNA_boxplot.pdf');

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
saveas(gca,'figures/benchmarkOriginalVsNewFPA_boxplot.pdf');
%% box plot of individual tissues 
% brain is predicted by rxn FPA (i = 3) (other examples: pancreas i = 12; Adrenal Gland i = 1,skin i =5; liver i = 6, muscle i = 10, Thyroid i = 16)
i = 16;

measuredTissue(i)
sum(logical(refMat(:,i)))

p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP2 = predMat2(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP2;predMat_all2(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched2'},length(enrichedMet_rFP2),1);repmat({'all2'},size(predMat_all2,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{met}';'all metabolites^{met}'}];
p(1) = ranksum(enrichedMet_rFP2,predMat_all2(:,i),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP;predMat_all(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{met+rxn}';'all metabolites^{met+rxn}'}];
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
saveas(gca,['figures/benchmarkMetabolomicsPredictor_boxplot_individualTissue_',measuredTissue{i},'.pdf']);

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

enrichedMet_rFP_noNetwork = predMat_noNetwork(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP_noNetwork;predMat_all_noNetwork(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_noNetwork'},length(enrichedMet_rFP_noNetwork),1);repmat({'all_noNetwork'},size(predMat_all_noNetwork,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{tranporter+rxn expression}';'all metabolites^{tranporter+rxn expression}'}];
p(2) = ranksum(enrichedMet_rFP_noNetwork,predMat_all_noNetwork(:,i),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat(:,i)),i);
boxData = [boxData;enrichedMet_rFP;predMat_all(:,i)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1),1)]];
xLabels = [xLabels;{'enriched metabolites^{FPA}';'all metabolites^{FPA}'}];
p(3) = ranksum(enrichedMet_rFP,predMat_all(:,i),  'Tail','right');


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
saveas(gca,['figures/benchmarkIntegrationBenefit_boxplot_individualTissue_',measuredTissue{i},'.pdf']);

%% box plot of individual tissues 
% i = 1, Adrenal Gland' RNA outperforms protein; also i = 2, Artery (but
% only 2 mets)
% i = 5, skin ; protein is better; i = 8 intestine 
i = 16;

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
saveas(gca,['figures/benchmarkProteinVsRNA_boxplot_individualTissue_',measuredTissue{i},'.pdf']);