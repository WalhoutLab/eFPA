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
%% load datasets - CHECK ONE EXAMPLE
nameMat = cell(size(refMat,1),size(refMat,2));
for i = 1:size(refMat,1)
    for j = 1:size(refMat,2)
        nameMat{i,j} = [allCmp_iHumanName{i},'-',measuredTissue{j}];
    end
end
MetLabels = nameMat(refMat == 1);
dorders = 0:20;
cutoffs = -1:0.01:1;
for dorder = dorders
    % protein prediction
    load(['output/FPA_HMDB_protein_TS_common_newFPA_weightedDist_order',num2str(dorder),'_tissueNetwork.mat']);

    % the internal rxn objectives 
    relFP = [output_collections.INT{1};output_collections.INT{2}];
    rowlabels = [cellfun(@(x) [x,'_f'],targetRxns_INT,'UniformOutput',false);
                  cellfun(@(x) [x,'_r'],targetRxns_INT,'UniformOutput',false)];   
    rmInd = all(isnan(relFP),2);
    relFP(rmInd,:) = [];
    rowlabels(rmInd) = [];
    relFP_INT = relFP;
    rowlabels_INT = rowlabels;

    % the transporter objectives
    relFP = [output_collections.TSP{1};output_collections.TSP{2}];
    rowlabels = [cellfun(@(x) [x,'_f'],targetRxns_TSP,'UniformOutput',false);
                  cellfun(@(x) [x,'_r'],targetRxns_TSP,'UniformOutput',false)];   
    rmInd = all(isnan(relFP),2);
    relFP(rmInd,:) = [];
    rowlabels(rmInd) = [];
    relFP_TSP = relFP;
    rowlabels_TSP = rowlabels;

    % the DM objectives
    relFP = [output_collections.DM{1};output_collections.DM{2}];
    rowlabels_DM = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
                  cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
    rmInd = all(isnan(relFP),2);
    relFP(rmInd,:) = [];
    rowlabels_DM(rmInd) = [];
    relFP_DM = relFP;
    %% make the nearest reaction potential 
    relFP_INT_ctd = normalize(relFP_INT,2,'center','median');
    %load('allCmp_iHumanName.mat');
    %allCmp_iHumanName = unique(allCmp_iHumanName);

    relFP_nearest = [];
    metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
    rowlabels_ID = regexprep(rowlabels_INT,'_.$','');
    S_logical = full(logical(model.S~=0));
    for i = 1:length(allCmp_iHumanName)
        myMet = allCmp_iHumanName{i};
        metInd = strcmp(myMet,metNames);
        myRxns = model.rxns(any(S_logical(metInd,:),1));
        myRxns = intersect(myRxns,targetRxns_INT);
        myInd = ismember(rowlabels_ID,myRxns);
        if any(myInd)
             relFP_nearest(i,:) = max(relFP_INT_ctd(myInd,:),[],1);
        else
            relFP_nearest(i,:)= zeros(1,size(relFP_nearest,2));
            fprintf('%s is not predicted (no internal rxn found)\n',allCmp_iHumanName{i});
        end
    end
    cmpName_nearest = allCmp_iHumanName;

    cmpName_nearest_network = cmpName_nearest;
    relFP_nearest_network = relFP_nearest;
    %% delta rFP matrix - transporter or DM or network-FPA-nearest - protein 
    relFP_sel = zeros(size(relFP_TSP,1),length(measuredTissue));
    relFP_TSP_ctd = normalize(relFP_TSP,2,'center','median');
    for i = 1:length(measuredTissue)
        relFP_sel(:,i) = max(relFP_TSP_ctd(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
    end
    relFP_sel_DM = zeros(size(relFP_DM,1),length(measuredTissue));
    relFP_wtd_ctd_DM = normalize(relFP_DM,2,'center','median');
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
        myInd = ismember(regexprep(rowlabels_TSP,'_.$',''),myRxns);

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
        myInd = ismember(regexprep(rowlabels_TSP,'_.$',''),myRxns);

        DMrxn = {['NewMet_',allCmp_iHumanName{j},'_f'],['NewMet_',allCmp_iHumanName{j},'_r']};
        DMInd = ismember(rowlabels_DM,DMrxn);

        if any(myInd) || any(DMInd)
            predMat2(j,:) = max([relFP_sel(myInd,:);relFP_sel_DM(DMInd,:)],[],1);
        else
            predMat2(j,:) = zeros(1,size(predMat2,2));
            fprintf('%s is not predicted\n',allCmp_iHumanName{j});
        end
    end

    %% we only look at recall, so we reshape the matrix to only look at positive tissue-metabolite pairs
    deltaFPmat(:,dorder+1) = predMat(refMat == 1);
    
    %% enrichment 
    for i = 1: length(cutoffs)
        predictions = predMat >= cutoffs(i);
        % test the narrow hypothesis 
        x = sum(predictions & refMat,1);
        M = size(predictions,1);
        K = sum(refMat,1);
        N = sum(predictions,1);
        p_enriMat(i,dorder+1) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    end
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

for i = 1: length(cutoffs)
    predictions = predMat >= cutoffs(i);
    % test the narrow hypothesis 
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_enri_base2(i,1) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
end
%%
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


deltaFPmat = [predMat_noNetwork(refMat == 1),predMat(refMat == 1), deltaFPmat];
IDs = [{'no integration','base 2 order 6'},strsplit(num2str(dorders))];

for i = 1: length(cutoffs)
    predictions = predMat_noNetwork >= cutoffs(i);
    % test the narrow hypothesis 
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_enri_noNet(i,1) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
end
%% save data
t = array2table(deltaFPmat);
t.Properties.RowNames = MetLabels;
writetable(t,'output/heatmapTbl.csv','WriteRowNames',1);
boundaries =  cell2table(IDs');
writetable(boundaries,'output/heatmapTbl_boundaries.csv');


p_enriMat = [p_enri_noNet,p_enri_base2, p_enriMat];
t = array2table(p_enriMat);
t.Properties.RowNames = strsplit(num2str(cutoffs));
writetable(t,'output/heatmapTbl_pEnri.csv','WriteRowNames',1);
boundaries =  cell2table(IDs');
writetable(boundaries,'output/heatmapTbl_boundaries_pEnri.csv');
%%
distMethod = 'euclidean';
%deltaFPmat_norm = deltaFPmat ./ max(abs(deltaFPmat),[],2); % relative r 
%deltaFPmat_norm(isnan(deltaFPmat_norm)) = 0;

cgo=clustergram(deltaFPmat(:,1:end),'RowLabels',MetLabels,'ColumnLabels',IDs(1:end),'RowPDist',distMethod,'Cluster', 'Column');
c=get(cgo,'ColorMap');
n = 100;
tmp = [zeros(n,1), linspace(1,0,n)',zeros(n,1)];
tmp = tmp(2:end,:);
cpr=[tmp; linspace(0,1,n)',zeros(n,1),zeros(n,1)
    ];
% cpr = cpr(199:-1:1,:);
set(cgo,'ColorMap',cpr);
set(cgo,'Symmetric',true);

set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 12)
%%
distMethod = 'euclidean';
p_enriMat_log = -log10(p_enriMat);

cgo=HeatMap(p_enriMat_log(:,1:end),'ColumnLabels',IDs(1:end),'RowLabels',strsplit(num2str(cutoffs)),'Symmetric','false');
    

