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
%% delta rFP matrix - transporter or DM or network-FPA-nearest - protein 
relFP_sel = zeros(size(relFP_wtd,1),length(measuredTissue));
relFP_wtd_ctd = normalize(relFP_wtd,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel(:,i) = max(relFP_wtd_ctd(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_DM = zeros(size(relFP_wtd_DM,1),length(measuredTissue));
relFP_wtd_ctd_DM = normalize(relFP_wtd_DM,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_DM(:,i) = max(relFP_wtd_ctd_DM(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_nearest = zeros(size(relFP_nearest_network,1),length(measuredTissue));
relFP_wtd_ctd_nearest = relFP_nearest_network;
for i = 1:length(measuredTissue)
    relFP_sel_nearest(:,i) = max(relFP_wtd_ctd_nearest(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

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

% delta rFP matrix - transporter or DM - protein 
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
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - transporter or DM or network-FPA-nearest - protein - original MERGE
relFP_sel_oriMERGE = zeros(size(relFP_wtd_oriMERGE,1),length(measuredTissue));
relFP_wtd_ctd_oriMERGE = normalize(relFP_wtd_oriMERGE,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_oriMERGE(:,i) = max(relFP_wtd_ctd_oriMERGE(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_DM_oriMERGE = zeros(size(relFP_wtd_DM_oriMERGE,1),length(measuredTissue));
relFP_wtd_ctd_DM_oriMERGE = normalize(relFP_wtd_DM_oriMERGE,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_DM_oriMERGE(:,i) = max(relFP_wtd_ctd_DM_oriMERGE(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
relFP_sel_nearest_oriMERGE = zeros(size(relFP_nearest_network_oriMERGE,1),length(measuredTissue));
relFP_wtd_ctd_nearest_oriMERGE = relFP_nearest_network_oriMERGE;
for i = 1:length(measuredTissue)
    relFP_sel_nearest_oriMERGE(:,i) = max(relFP_wtd_ctd_nearest_oriMERGE(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

predMat_oriMERGE = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_oriMERGE,'_.$',''),myRxns);
    
    DMrxn = {['NewMet_',allCmp_iHumanName{j},'_f'],['NewMet_',allCmp_iHumanName{j},'_r']};
    DMInd = ismember(rowlabels_DM_oriMERGE,DMrxn);
    
    NearestInd = ismember(cmpName_nearest_network_oriMERGE,allCmp_iHumanName{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_oriMERGE(j,:) = max([relFP_sel_oriMERGE(myInd,:);relFP_sel_DM_oriMERGE(DMInd,:);relFP_sel_nearest_oriMERGE(NearestInd,:)],[],1);
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

%% delta rFP matrix - transporter or DM or network-FPA-nearest - RNA all
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

%% delta rFP matrix - transporter or DM - no network integration
relFP_sel_noNetwork = zeros(size(relFP_wtd_noNetwork,1),length(measuredTissue));
relFP_wtd_ctd_noNetwork = normalize(relFP_wtd_noNetwork,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_noNetwork(:,i) = max(relFP_wtd_ctd_noNetwork(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

relFP_sel_DM_noNetwork = zeros(size(relFP_nearest,1),length(measuredTissue));
relFP_wtd_ctd_DM_NoNetwork = relFP_nearest;
for i = 1:length(measuredTissue)
    relFP_sel_DM_noNetwork(:,i) = max(relFP_wtd_ctd_DM_NoNetwork(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end

predMat_noNetwork = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_noNetwork,'_.$',''),myRxns);
    
    DMInd = ismember(cmpName_nearest,allCmp_iHumanName{j});
    if any(myInd) || any(DMInd)
        predMat_noNetwork(j,:) = max([relFP_sel_noNetwork(myInd,:);relFP_sel_DM_noNetwork(DMInd,:)],[],1);
    else
        predMat_noNetwork(j,:) = zeros(1,size(predMat_noNetwork,2));
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
%% GM1 and thiamine
% brain Ganglioside Gm1 GA2 digalactosylceramide
% Muscle thiamin
metname = 'Ganglioside Gm1';
tissuename = 'Brain';
ind1 = strcmp(allCmp_iHumanName,metname);
ind2 = strcmp(measuredTissue,tissuename);
refMat(ind1,ind2)
proteinPred = predMat(ind1,ind2)
RNAPred = predMat_RNAcomm(ind1,ind2)
RNAPred_all = predMat_RNAall(ind1,ind2)
figure(1)
c = categorical({'Protein','RNA','RNA (all genes)'});
bar(c, [proteinPred,RNAPred,RNAPred_all],'FaceColor','#808080')
ylabel(['delta rFP '])
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Interpreter = 'none';
plt.TickDir = 'out';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.Title = [metname,' - ',tissuename];
plt.export(['figures/HMDB_examples_',metname,'_',tissuename,'.pdf']);
