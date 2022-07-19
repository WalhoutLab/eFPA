%% About 
% systematically benchmark FPA (the principle of local expression
% dictating metabolic flux/functions) and compare protein data with RNA
% data as the input, based on HMDB reference tissue-enriched metabolite
% dataset
%% set up environments
setEnvForAnalysis
addpath('./PlotPub/lib')

% define reaction sets
% the exchange reactions with environment 
excRxns = model.rxns(findExcRxns_XL(model));
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
excRxns = intersect(excRxns,EXinvolvedRxns);

% the transporter (with env) reactions
allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
TSP = [];
for i = 1:length(allCmp_iHumanName)
    myMet_e = {[allCmp_iHumanName{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,allCmp_iHumanName(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    % we define the transporter as the reactions that contain the same
    % metabolite in [e] and another compartment (cellular) in the same
    % reaction
    candidate = intersect(myRxns_non_e,myRxns_e);
    % check if is on diff side of the reaction
    if ~isempty(candidate)
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                TSP = union(TSP,candidate(j));
            end
        end
    end
end
% some special transporter will be missed, we add back 
envTspRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
envTspRxns = intersect(envTspRxns,EXinvolvedRxns);
tspRxns = union(TSP, envTspRxns);

% the internal transporters 
% internal transporters are hard to define, since a reaction can span two
% compartments internally but it is not a real transporter. So, we use the
% subsys annotation as a compromise
intTspRxns = setdiff(model.rxns(ismember(model.subSystems,{'Transport reactions'})),tspRxns);

% internal regular reactions
intRxns = setdiff(model.rxns, [excRxns; tspRxns; intTspRxns]);

% regular met-analysis target (transporters)
targetRxns = tspRxns; % in this script, the targetrxns by default refers to the transporter objectives

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
cmpName_nearest_noNetwork = cmpName_nearest;
relFP_nearest_noNetwork = relFP_nearest;
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
cutoffs = -1.01:0.01:1.01;

% protein - common genes - this is considered as default 
OverallPrecision =[];
OverallRecall =[];
TPR = [];
FPR = [];
% protein - common genes - without asso. int. rxns
OverallPrecision2 =[];
OverallRecall2 =[];
TPR2 = [];
FPR2 = [];

OverallPrecision_RNAcomm =[];
OverallRecall_RNAcomm =[];
TPR_RNAcomm = [];
FPR_RNAcomm = [];

OverallPrecision_RNAall =[];
OverallRecall_RNAall =[];
TPR_RNAall = [];
FPR_RNAall = [];

OverallPrecision_noNetwork =[];
OverallRecall_noNetwork =[];
TPR_noNetwork = [];
FPR_noNetwork = [];

OverallPrecision_noNetwork_tsp =[];
OverallRecall_noNetwork_tsp =[];
TPR_noNetwork_tsp = [];
FPR_noNetwork_tsp = [];

for i = 1: length(cutoffs)
    predictions = predMat >= cutoffs(i) & predMat < 1.01;% if it is larger than 1, we basically consider it numerical error
    % although the numerical-error related predictions are likely to be a
    % real high delta rFP, we exclude it to be conservative. Since human1
    % model is too large, few numeric errors in a large-scale computation
    % is hard to avoid. 
    
    % PR curve
    TP = sum(predictions & refMat,1);
    FP = sum(predictions & ~refMat,1);
    FN = sum(~predictions & refMat,1);
    TN = sum(~predictions & ~refMat,1);
    OverallPrecision(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall(i) = sum(TP) ./ (sum(TP) + sum(FN));
    TPR(i) = sum(TP) ./ (sum(TP) + sum(FN));
    FPR(i) = sum(FP) ./ (sum(TN) + sum(FP));
    
    
    % apply the same calculation to other predictions
    predictions = predMat2 >= cutoffs(i) & predMat2 < 1.01;
    TP = sum(predictions & refMat,1);
    FP = sum(predictions & ~refMat,1);
    FN = sum(~predictions & refMat,1);
    TN = sum(~predictions & ~refMat,1);
    OverallPrecision2(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall2(i) = sum(TP) ./ (sum(TP) + sum(FN));
    TPR2(i) = sum(TP) ./ (sum(TP) + sum(FN));
    FPR2(i) = sum(FP) ./ (sum(TN) + sum(FP));
    
    predictions = predMat_RNAcomm >= cutoffs(i) & predMat_RNAcomm < 1.01;
    TP = sum(predictions & refMat,1);
    FP = sum(predictions & ~refMat,1);
    FN = sum(~predictions & refMat,1);
    TN = sum(~predictions & ~refMat,1);
    OverallPrecision_RNAcomm(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall_RNAcomm(i) = sum(TP) ./ (sum(TP) + sum(FN));
    TPR_RNAcomm(i) = sum(TP) ./ (sum(TP) + sum(FN));
    FPR_RNAcomm(i) = sum(FP) ./ (sum(TN) + sum(FP));
    
   
    predictions = predMat_RNAall >= cutoffs(i) & predMat_RNAall < 1.01;
    TP = sum(predictions & refMat,1);
    FP = sum(predictions & ~refMat,1);
    FN = sum(~predictions & refMat,1);
    TN = sum(~predictions & ~refMat,1);
    OverallPrecision_RNAall(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall_RNAall(i) = sum(TP) ./ (sum(TP) + sum(FN));
    TPR_RNAall(i) = sum(TP) ./ (sum(TP) + sum(FN));
    FPR_RNAall(i) = sum(FP) ./ (sum(TN) + sum(FP));

    
    predictions = predMat_noNetwork >= cutoffs(i) & predMat_noNetwork < 1.01;
    TP = sum(predictions & refMat,1);
    FP = sum(predictions & ~refMat,1);
    FN = sum(~predictions & refMat,1);
    TN = sum(~predictions & ~refMat,1);
    OverallPrecision_noNetwork(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall_noNetwork(i) = sum(TP) ./ (sum(TP) + sum(FN));
    TPR_noNetwork(i) = sum(TP) ./ (sum(TP) + sum(FN));
    FPR_noNetwork(i) = sum(FP) ./ (sum(TN) + sum(FP));

    
    predictions = predMat_noNetwork_tsp >= cutoffs(i) & predMat_noNetwork_tsp < 1.01;
    TP = sum(predictions & refMat,1);
    FP = sum(predictions & ~refMat,1);
    FN = sum(~predictions & refMat,1);
    TN = sum(~predictions & ~refMat,1);
    OverallPrecision_noNetwork_tsp(i) = sum(TP) ./ (sum(TP) + sum(FP));
    OverallRecall_noNetwork_tsp(i) = sum(TP) ./ (sum(TP) + sum(FN));
    TPR_noNetwork_tsp(i) = sum(TP) ./ (sum(TP) + sum(FN));
    FPR_noNetwork_tsp(i) = sum(FP) ./ (sum(TN) + sum(FP));
    
end
OverallPrecision(isnan(OverallPrecision)) =1;
OverallPrecision2(isnan(OverallPrecision2)) =1;
OverallPrecision_RNAall(isnan(OverallPrecision_RNAall)) =1;
OverallPrecision_noNetwork(isnan(OverallPrecision_noNetwork)) =1;
OverallPrecision_noNetwork_tsp(isnan(OverallPrecision_noNetwork_tsp)) =1;

figure(1)
hold on
plot(OverallRecall, OverallPrecision,'.-');
hold off
ylabel('Precison')
xlabel('recall')

figure(2)
hold on
plot(FPR, TPR,'.-');
plot([0 1],[0,1],'-')
hold off
ylabel('TPR')
xlabel('FPR')

figure(3)
hold on
plot(OverallRecall2, OverallPrecision2,'-');
plot(OverallRecall, OverallPrecision,'-');
yline(sum(sum(refMat))./(size(refMat,1).*size(refMat,2)),'--','LineWidth',2,'Color',[0.5 0.5 0.5])
hold off
legend({'FPA - exchaning','FPA - exchaning+internal'});
xlabel('Recall')
ylabel('Precision')
ylim([0.14,0.3]);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3, 3];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.TickDir = 'out';
plt.AxisLineWidth = 1;
plt.export('figures/NoTrack_benchmarkMetabolomicsPredictor_PRC.pdf');

figure(4)
hold on
plot(OverallRecall, OverallPrecision,'-');
plot(OverallRecall_noNetwork_tsp, OverallPrecision_noNetwork_tsp,'-');
plot(OverallRecall_noNetwork, OverallPrecision_noNetwork,'-');
yline(sum(sum(refMat))./(size(refMat,1).*size(refMat,2)),'--','LineWidth',2,'Color',[0.5 0.5 0.5])
hold off
ylim([0.14,0.3]);
legend({'FPA','Transporter expression',sprintf('Expression of transporter \nand associated rxn(s)')});
xlabel('Recall')
ylabel('Precision')
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [3, 3];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.TickDir = 'out';
plt.AxisLineWidth = 1;
plt.export('figures/NoTrack_benchmarkIntegrationBenefit_PRC.pdf');

figure(5)
[A B] = sort(OverallRecall);
AUC1 = trapz(OverallRecall(B), OverallPrecision(B));
[A B] = sort(OverallRecall_noNetwork_tsp);
AUC2 = trapz(OverallRecall_noNetwork_tsp(B), OverallPrecision_noNetwork_tsp(B));
[A B] = sort(OverallRecall_noNetwork);
AUC3 = trapz(OverallRecall_noNetwork(B), OverallPrecision_noNetwork(B));
c = categorical({'FPA','Transporter expression',sprintf('Expression of transporter \nand associated rxn(s)')});
c = reordercats(c,{'Transporter expression',sprintf('Expression of transporter \nand associated rxn(s)'),'FPA'});
hold on 
bar(c, [AUC1,AUC2,AUC3],'FaceColor','#808080')
ylabel('AUPRC')
yline(sum(sum(refMat))./(size(refMat,1).*size(refMat,2)),'--','LineWidth',2,'Color',[0.5 0.5 0.5])
hold off
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [3, 3];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.TickDir = 'out';
plt.AxisLineWidth = 1;
plt.export('figures/NoTrack_benchmarkIntegrationBenefit_AUPRC.pdf');


figure(6)
hold on
plot(OverallRecall, OverallPrecision,'-');
plot(OverallRecall, OverallPrecision_RNAcomm,'-');
plot(OverallRecall, OverallPrecision_RNAall,'-');
yline(sum(sum(refMat))./(size(refMat,1).*size(refMat,2)),'--','LineWidth',2,'Color',[0.5 0.5 0.5])
hold off
ylim([0.14,0.3]);
legend({'Protein FPA','RNA FPA','RNA FPA (all genes)'});
xlabel('Recall')
ylabel('Precision')
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [3, 3];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.AxisLineWidth = 1;
plt.TickDir = 'out';
%plt.export('figures/NoTrack_benchmarkProteinVsRNA.pdf');


figure(7)
hold on
plot(FPR, TPR,'-');
plot(FPR_noNetwork_tsp, TPR_noNetwork_tsp,'-');
plot(FPR_noNetwork, TPR_noNetwork,'-');
hold off
legend({'FPA','Transporter expression',sprintf('Expression of transporter \nand associated rxn(s)')});
xlabel('FPR')
ylabel('TPR')
hline = refline([1 0]);
hline.Color = 'k';
hline.LineStyle = ':';
hline.HandleVisibility = 'off';
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [3, 3];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.TickDir = 'out';
plt.AxisLineWidth = 1;
plt.export('figures/NoTrack_benchmarkIntegrationBenefit_ROC.pdf');

