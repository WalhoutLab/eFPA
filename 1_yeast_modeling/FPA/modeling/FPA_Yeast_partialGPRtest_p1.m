%% notes
% flux is in  (moles / hr / mL cells); could be further normalized to
% mmole/hr/gDW by chemostat info: gDCW/ml. Should try both and see
% correlation change
addpath /share/pkg/gurobi/900/linux64/matlab/
addpath ~/cobratoolbox/
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/
initCobraToolbox(false)

%% 1. load the model and prepare the model
addpath('./../scripts/')
model = loadYeatModel();
% the following nutrients need to be set manually
% to start with basic FPA, we allow unlimited exchange
% phosphate exchange
model.lb(strcmp(model.rxns,'r_2005')) = -1000;
% glucose exchange
model.lb(strcmp(model.rxns,'r_1714')) = -1000;
% ammonium exchange 
model.lb(strcmp(model.rxns,'r_1654')) = -1000;
% uracil 
model.lb(strcmp(model.rxns,'r_2090')) = -1000;
% leucine
model.lb(strcmp(model.rxns,'r_1899')) = -1000;
% maintanence 
model = changeRxnBounds(model,'r_4046',0,'l'); % maintance 
model = changeRxnBounds(model,'r_4046',1000,'u'); % maintance 

% allow all exchange
% model.lb(findExcRxns(model)) = -1000;

%% 3. load the expression files, distance matrix, and other optional inputs
% load the distance matrix
% users can uncomment the following codes to load from distance calculator
% output; here we load directly from saved matlab variable because of file
% size restriction of GitHub
distance_raw = readtable('./../input/YeastJoshua/distanceMatrix_weighted.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github
labels = distance_raw.Properties.VariableNames;
labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
distMat_raw = table2array(distance_raw);
distMat_min = zeros(size(distance_raw,1),size(distance_raw,2));
for i = 1:size(distMat_min,1)
    for j = 1:size(distMat_min,2)
        distMat_min(i,j) = min([distMat_raw(i,j),distMat_raw(j,i)]);
    end
end
distMat_wtd = distMat_min;

distance_raw = readtable('./../input/YeastJoshua/distanceMatrix.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github
labels = distance_raw.Properties.VariableNames;
labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
distMat_raw = table2array(distance_raw);
distMat_min = zeros(size(distance_raw,1),size(distance_raw,2));
for i = 1:size(distMat_min,1)
    for j = 1:size(distMat_min,2)
        distMat_min(i,j) = min([distMat_raw(i,j),distMat_raw(j,i)]);
    end
end
distMat_normal = distMat_min;

% load the special penalties 
% In general, we recommend tp set penalty for all Exchange, Demand,
% and Sink reactions to 0 to not penaltize the external reactions. Users 
% may need to interactively tune their special penalties for best flux
% distribution in the FPA calculation
extRxns = model.rxns(findExcRxns(model));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

% we dont recomand any specific special distance for generic model; In the
% dual model, the special distance was used to discourage using of side
% metabolites. Since side/storage metabolites are not applicable for
% generic model, we don't use any special distance. 
% manualDist = {};

proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

% make the blocklist to represent the condition-specific nutrient 
p_lim = {'r_1899_r','r_2090_r'};
c_lim = {'r_1899_r','r_2090_r'};
n_lim = {'r_1899_r','r_2090_r'};
l_lim = {'r_2090_r'};
u_lim = {'r_1899_r'};
blocklist = [repmat({p_lim},1,5),repmat({c_lim},1,5),repmat({n_lim},1,5),repmat({l_lim},1,5),repmat({u_lim},1,5),{{}}];
%% flux table 
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;
% normalize flux
% flux is in  (moles / hr / mL cells); could be further normalized to
% mmole/hr/gDW by chemostat info: gDCW/ml. Should try both and see
% correlation change
%fluxMat_normalized = fluxMat;

% dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
% factor = repmat(dwTbl.gDCW_mL',size(fluxMat,1),1);
% fluxMat_normalized = fluxMat_normalized * 1000 ./ factor; %mmoles/hr/gDW

%GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
%fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
%  dilutionFactor = repmat([0.05 0.1 0.16 0.22 0.30],size(fluxMat,1),5);
%  fluxMat_double_normalized = fluxMat_normalized ./ dilutionFactor;
%% load merged levels
load('output/normalizedLevels.mat');
%% precalc penalty 
penalty_pro = ones(length(model.rxns),size(fluxMat,2)+1);
[A B] = ismember(model.rxns,valid_rxns_pro_perPro);
penalty_pro(A,1:(end-1)) = ones(size(normalizedLevel_pro_perPro(B(A),:),1),size(normalizedLevel_pro_perPro(B(A),:),2)) ./ normalizedLevel_pro_perPro(B(A),:);

penalty_pro_raw = ones(length(model.rxns),size(fluxMat,2)+1);
[A B] = ismember(model.rxns,valid_rxns_pro_perDW);
penalty_pro_raw(A,1:(end-1)) = ones(size(normalizedLevel_pro_perDW(B(A),:),1),size(normalizedLevel_pro_perDW(B(A),:),2)) ./ normalizedLevel_pro_perDW(B(A),:);

% apply additional penalty to the exchange of limiting nutrients
for i = 1:size(manualPenalty,1)
    penalty_pro(strcmp(model.rxns,manualPenalty{i,1}),:) = manualPenalty{i,2};
end
penalty_pro(strcmp(model.rxns,'r_2005'),1:5) = 10;
penalty_pro(strcmp(model.rxns,'r_1714'),6:10) = 10;
penalty_pro(strcmp(model.rxns,'r_1654'),11:15) = 10;

for i = 1:size(manualPenalty,1)
    penalty_pro_raw(strcmp(model.rxns,manualPenalty{i,1}),:) = manualPenalty{i,2};
end
penalty_pro_raw(strcmp(model.rxns,'r_2005'),1:5) = 10;
penalty_pro_raw(strcmp(model.rxns,'r_1714'),6:10) = 10;
penalty_pro_raw(strcmp(model.rxns,'r_1654'),11:15) = 10;



%%
myCluster = parcluster('local');
myCluster.NumWorkers = 128;
saveProfile(myCluster);
myCluster
parpool(39,'SpmdEnabled',false);
%% original FPA with original distance - normalized ones
% setup some basic parameters for FPA
n = [0:0.5:5 6 7 8 9 10 15 20]; % note: each will take 10GB mem! dont titrate all "n" and target rxns together!
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
targetRxns = rxnLabel;
% targetRxns = {'r_4041'};
%% FPA2 (new decay formula): orginal distance - normalized ones
% exponential filtering could strongly filter out the influence from distal
% distance expression; and it could achive more binary filtering (lower
% than cutoff ==> strongly used; higher than distance cutoff ==> strongly
% filtered)

% by tuning the base of the exponential formula, we can define how binary
% the transition should be 
n2 = [0 1 2 2.5 3 3.5 4 5 6 7 8 9 10 20 50];
base = [2 100];
%% FPA2 (new decay formula): wtd distance - normalized ones
[FP_collection_2] = FPA2_cluster_multiPara(model,targetRxns,{},distMat_wtd,labels,n2, manualPenalty,{},max(distMat_wtd(~isinf(distMat_wtd))),blocklist, {},true,penalty_pro,base);
save('output/partialGPR/Titration_relativeExp_wtdDist_expDecay_NoExcluded.mat','FP_collection_2','targetRxns','n2','base');

