%% About
% run the systematic FPA analysis on 232 yeast reactions. this is a large
% scale computation that requires a lab server to run
% IMPORTANT NOTICE:
% Please move this script to the /1_yeast_modeling/FPA/modeling/ before
% running! Otherwise there would be path issues!
%% paths
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
%% 3. load the expression files, distance matrix, and other optional inputs
% load the distance matrix
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
extRxns = model.rxns(findExcRxns(model));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

% load condition names
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
%% setup some basic parameters for FPA
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
targetRxns = rxnLabel;
n2 = 0:0.5:40;
alpha = 1;
nsample = 1000;
seed = 1126;
rng(seed);
mkdir(['output/randomization_simpleDecay/seed_',num2str(seed),'/'])
%% complete randomization
for c = 1:(nsample+1)
    %% load expression levels
    load('output/normalizedLevels_partialExcluded.mat');
    %% shuffling the gene labels 
    if c == 1
        % dont shuffle as control
    else
        valid_rxns_pro_perPro = valid_rxns_pro_perPro(randperm(length(valid_rxns_pro_perPro)));
    end
    %% precalc penalty 
    penalty_pro = ones(length(model.rxns),size(fluxMat,2)+1);
    [A B] = ismember(model.rxns,valid_rxns_pro_perPro);
    penalty_pro(A,1:(end-1)) = ones(size(normalizedLevel_pro_perPro(B(A),:),1),size(normalizedLevel_pro_perPro(B(A),:),2)) ./ normalizedLevel_pro_perPro(B(A),:);

    % apply additional penalty to the exchange of limiting nutrients
    for i = 1:size(manualPenalty,1)
        penalty_pro(strcmp(model.rxns,manualPenalty{i,1}),:) = manualPenalty{i,2};
    end
    penalty_pro(strcmp(model.rxns,'r_2005'),1:5) = 10;
    penalty_pro(strcmp(model.rxns,'r_1714'),6:10) = 10;
    penalty_pro(strcmp(model.rxns,'r_1654'),11:15) = 10;
    %% FPA 
    [FP_collection_2] = FPA2_clusterWrapper_multiPara(model,targetRxns,{},distMat_wtd,labels,6, manualPenalty,{},max(distMat_wtd(~isinf(distMat_wtd))),blocklist, {},false,penalty_pro,alpha,2);
    save(['output/randomization_simpleDecay/seed_',num2str(seed),'/rand_',num2str(c),'_default_FPA.mat'],'FP_collection_2','targetRxns');

    [FP_collection_2] = FPA2_clusterWrapper_simpleDecay(model,targetRxns,{},distMat_wtd,labels,n2, manualPenalty,{},max(distMat_wtd(~isinf(distMat_wtd))),blocklist, {},false,penalty_pro,alpha);
    save(['output/randomization_simpleDecay/seed_',num2str(seed),'/rand_',num2str(c),'_flexi_FPA.mat'],'FP_collection_2','targetRxns');
end

