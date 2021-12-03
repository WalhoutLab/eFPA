%% About
% run yeast FPA analysis locally on the laptop for small scale analysis 
%%
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/
initCobraToolbox
%% 1. load the model and prepare the model
addpath('./../scripts/')
model = loadYeatModel();
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
distMat = distMat_min;
% set the special penalties 
% In general, we recommend tp set penalty for all Exchange, Demand,
% and Sink reactions to 0 to not penaltize the external reactions. Users 
% may need to interactively tune their special penalties for best flux
% distribution in the FPA calculation
extRxns = model.rxns(findExcRxns(model));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

% load protein levels (this is just to obtain the condition names)
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

% make the blocklist to represent the condition-specific nutrient 
p_lim = {'r_1899_r','r_2090_r'};
c_lim = {'r_1899_r','r_2090_r'};
n_lim = {'r_1899_r','r_2090_r'};
l_lim = {'r_2090_r'};
u_lim = {'r_1899_r'};
blocklist = [repmat({p_lim},1,5),repmat({c_lim},1,5),repmat({n_lim},1,5),repmat({l_lim},1,5),repmat({u_lim},1,5),{{}}];
%% load flux table 
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;
% normalize flux
fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
labels2 = regexprep(conditions,'_','-');
%% load expression levels
load('output/normalizedLevels_partialExcluded.mat');
%% precalc penalty (1/normalized expression coefficient)
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
%% run the improved FPA
parpool(2)
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
%% run on an example reaction (large scale running is done in other scripts)
targetRxn = {'r_0565'};
distBound = 6;
base = 2;
[FP_collection,FP_solutions_collection] = FPA2(model,targetRxn,{},distMat,labels,distBound, manualPenalty,{},max(distMat(~isinf(distMat))),blocklist, {},true,penalty_pro,base);
% calculate relative flux potential rFP
FP = FP_collection;
relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
i = find(strcmp(rxnLabel,targetRxn));
for j = 1:(size(FP,2)-1)
    if  mean(fluxMat(i,:)) > 0 
        if ~isnan(FP{1,j}(1))
            relFP(1,j) = FP{1,j}(1) ./ FP{1,end}(1);
        elseif ~isnan(FP{1,j}(2))
            relFP(1,j) = FP{1,j}(2) ./ FP{1,end}(2);
        else
            error('check');
        end
    else
        if ~isnan(FP{1,j}(2))
            relFP(1,j) = FP{1,j}(2) ./ FP{1,end}(2);
        elseif ~isnan(FP{1,j}(1))
            relFP(1,j) = FP{1,j}(1) ./ FP{1,end}(1);
        else
            error('check');
        end
    end
end
%% visualize correlation 
myLevel = relFP;
myFluxLevel = abs(fluxMat_normalized(strcmp(rxnLabel,targetRxn),:));

figure(1);
lm = fitlm(myLevel,myFluxLevel);
plot(lm);
text(myLevel,myFluxLevel,labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')
[corrR,corrP] = corr(myFluxLevel',myLevel')
%% check flux distribnution in the FPA prediction 
sample = 'C0_05';
myrxn = targetRxn;
if mean(fluxMat(i,:)) > 0 
    solution = FP_solutions_collection{strcmp(myrxn,targetRxn),strcmp(sample,conditions)}{1};
else
    solution = FP_solutions_collection{strcmp(myrxn,targetRxn),strcmp(sample,conditions)}{2};
end
mytbl = listFPtbl(model,solution);
%% track flux distribution by metabolite
model_irre = convertToIrreversible(model);
mytbl = listRxn(model_irre, solution.full,'s_1556[c]');
