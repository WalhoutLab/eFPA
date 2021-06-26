%% the p-value of boxplots should be changed to rank-sum test
%%
setEnvForAnalysis
addpath('PlotPub/lib')
load('HMDB_objectives.mat');
%% the HMDB reference dataset
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
allCmp_iHumanName_valid = allCmp_iHumanName;

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
%% load control - delta rFP matrix - transporter or DM - no network integration - and base2
% control - no network
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

allCmp_iHumanName2 = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName2 = setdiff(allCmp_iHumanName2,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
targetRxns_allMetDM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName2,'UniformOutput',false);

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

load output/FPA_transporter_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_noNetwork = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_noNetwork(rmInd) = [];
relFP_wtd_noNetwork = relFP;
load output/pseudoDM_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat

%
allCmp_iHumanName = allCmp_iHumanName_valid;
% delta rFP matrix - transporter or DM or network-FPA-nearest - protein 
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
predMat_detail = {};
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels,'_.$',''),myRxns);
    
    DMrxn = {['NewMet_',allCmp_iHumanName{j},'_f'],['NewMet_',allCmp_iHumanName{j},'_r']};
    DMInd = ismember(rowlabels_DM,DMrxn);
    
    NearestInd = ismember(cmpName_nearest_network,allCmp_iHumanName{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat(j,:) = max([relFP_sel(myInd,:);relFP_sel_DM(DMInd,:);relFP_sel_nearest(NearestInd,:)],[],1);
        detail = array2table([relFP_sel(myInd,:);relFP_sel_DM(DMInd,:);relFP_sel_nearest(NearestInd,:)]);
        detail.Properties.RowNames = [rowlabels(myInd);rowlabels_DM(DMInd);cmpName_nearest_network(NearestInd)];
        detail.Properties.VariableNames = measuredTissue;
        for i = 1:length(measuredTissue)
            predMat_detail(j,i) = {detail(:,i)};
        end
    else
        predMat(j,:) = zeros(1,size(predMat,2));
        predMat_detail(j,:) = {};
        fprintf('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

nameMat = cell(size(refMat,1),size(refMat,2));
for i = 1:size(refMat,1)
    for j = 1:size(refMat,2)
        nameMat{i,j} = [allCmp_iHumanName{i},'-',measuredTissue{j}];
    end
end
MetLabels = nameMat(refMat == 1);
deltaFPmat = predMat(refMat == 1);
deltaFPmat_detail = predMat_detail(refMat == 1);

%% save data
t = array2table(deltaFPmat);
t.Properties.RowNames = MetLabels;
t.detail = deltaFPmat_detail;
t.Properties.VariableNames = {'max_delta_rFP','delta_rFP_details'};
save('output/tissue_enriched_met_delta_rFP_tbl.mat','t')

% merge the same metabolite
allCmp_iHumanName_unique = unique(allCmp_iHumanName);
for i = 1:length(allCmp_iHumanName_unique)
    if sum(strcmp(allCmp_iHumanName,allCmp_iHumanName_unique{i})) ==1 
        predMat_uni(i,:) = predMat(strcmp(allCmp_iHumanName,allCmp_iHumanName_unique{i}),:);
    else
        if size(unique(predMat(strcmp(allCmp_iHumanName,allCmp_iHumanName_unique{i}),:),'rows'),1) > 1
            error('??');
        end
        predMat_uni(i,:) = unique(predMat(strcmp(allCmp_iHumanName,allCmp_iHumanName_unique{i}),:),'rows');
    end
end
delta_rFP_Mat = array2table(predMat_uni);
delta_rFP_Mat.Properties.RowNames = allCmp_iHumanName_unique;
delta_rFP_Mat.Properties.VariableNames = measuredTissue;
save('output/delta_rFP_tbl_HMDBmet.mat','delta_rFP_Mat')
