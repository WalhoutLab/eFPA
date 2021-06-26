%% first save rxn-centered tbl
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

load output/FPA_rxn_RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
relFP(rmInd,:) = [];
relFP_RNA = relFP;
%% thresholding
% cutoffs = 0:0.01:1;
% for i = 1: length(cutoffs)
%     N_pass_rna(i) = sum(sum(abs(relFP_RNA - median(relFP_RNA,2)) > cutoffs(i)));
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
relFP_filtered_RNA = relFP_RNA - median(relFP_RNA,2);

% keep = any(abs(relFP_filtered) > 0.2 | abs(relFP_filtered_RNA) > 0.2,2);

% relFP_filtered = relFP_filtered(keep,:);
% rowlabels_filtered = rowlabels_filtered(keep);
% relFP_filtered_RNA = relFP_filtered_RNA(keep,:);

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

lib = 'RNA';
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
    highestSites = conditions((relFP_filtered_RNA(j,:) > 0.2));
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
IsInOFD_RNA = IsInOFD;

sum(IsInOFD_pro)
sum(IsInOFD_RNA)
% label out the in OFD
for i = 1:length(rowlabels_filtered)
    if IsInOFD_pro(i)
        rowlabels_filtered{i} = [rowlabels_filtered{i},'*'];
    end
end
% relFP_filtered = relFP_filtered(IsInOFD_pro | IsInOFD_RNA,:);
% relFP_filtered_RNA = relFP_filtered_RNA(IsInOFD_pro | IsInOFD_RNA,:);
% rowlabels_filtered = rowlabels_filtered(IsInOFD_pro | IsInOFD_RNA,:);

%% save the table
IDs = regexprep(conditions,'_','-');
relFluxPotential = array2table(relFP_filtered);
relFluxPotential.Properties.VariableNames = IDs;
relFluxPotential.Properties.RowNames = rowlabels_filtered;
writetable(relFluxPotential,'output/supp_delta_rFP_rxn_protein.csv','WriteRowNames',true);

%% the metabolite table
clear
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
%% background delta rFP matrix - transporter or DM or nearest networkFPA - protein
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
    
    DMrxn = {['NewMet_',qryMets{j},'_f'],['NewMet_',qryMets{j},'_r']};
    DMInd = ismember(rowlabels_allMetDM,DMrxn);
    
    NearestInd = ismember(cmpName_nearest_network,qryMets{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_all = [predMat_all;max([relFP_sel(myInd,:);relFP_sel_allMetDM(DMInd,:);relFP_sel_nearest(NearestInd,:)],[],1)];
    else
        predMat_all = [predMat_all;nan(1,size(predMat_all,2))];
    end
end

predMat_all(abs(predMat_all)>1.01) = NaN;% 
%% save the table
relFluxPotential = array2table(predMat_all);
relFluxPotential.Properties.VariableNames = measuredTissue;
relFluxPotential.Properties.RowNames = qryMets;
writetable(relFluxPotential,'output/supp_delta_rFP_met_protein.csv','WriteRowNames',true);



