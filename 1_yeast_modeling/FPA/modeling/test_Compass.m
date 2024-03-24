function n_corr_COMPASS = test_Compass(n)
addpath ./../other_methods/COMPASS/
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
% load distance matrix - COMPASS DOES NOT NEED DISTANCE BUT WE USE AN FPA
% IMPLEMEMTATION SO JUST PASS THROUGH (THE DISTANCE ORDER OF ZERO WILL
% CANCEL OUT ALL DISTANCE EFFECTS)
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

% for the penalty matrix (and special constriants), we aligned them to eFPA
% so that they are more comparable. We can also remove these things later
% to see effects
extRxns = model.rxns(findExcRxns(model));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

% make the blocklist to represent the condition-specific nutrient 
p_lim = {'r_1899_r','r_2090_r'};
c_lim = {'r_1899_r','r_2090_r'};
n_lim = {'r_1899_r','r_2090_r'};
l_lim = {'r_2090_r'};
u_lim = {'r_1899_r'};
blocklist = [repmat({p_lim},1,5),repmat({c_lim},1,5),repmat({n_lim},1,5),repmat({l_lim},1,5),repmat({u_lim},1,5),{{}}];

fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;

% flux is in  (moles / hr / mL cells); no conversion is needed. 
% in fact, correlation got worse if we normalzie the flux to / gDW first!
fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);


% directly use the parsed expression levels from eFPA - this might be
% slightly different than compass but we have no way to use the original
% compass on a new model
load('output/normalizedLevels_partialExcluded.mat'); 
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

%% original FPA with original distance - normalized ones
% setup some basic parameters for FPA
% n = [0 2.5]; % n=0 is compass setup; 2.5 adds in a distance decay
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
targetRxns = rxnLabel;
% targetRxns = {'r_4041'};
%The FPA is designed with parfor loops for better speed, so we first initial the parpool
[FP_collection] = FPA_COMPASS_LIKE(model,targetRxns,{},distMat_normal,labels,n, manualPenalty,{},max(distMat_normal(~isinf(distMat_normal))),blocklist, {},true,penalty_pro);

for zz = 1:length(n)
    FP = FP_collection{zz};
    relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
    for i = 1:size(FP,1)
        for j = 1:(size(FP,2)-1)
            if  mean(fluxMat(i,:)) > 0 
                if ~isnan(FP{i,j}(1))
                    relFP(i,j) = -log(FP{i,j}(1)+1);
                elseif ~isnan(FP{i,j}(2))
                    relFP(i,j) = -log(FP{i,j}(2)+1);
                else
                    error('check');
                end
            else
                if ~isnan(FP{i,j}(2))
                    relFP(i,j) = -log(1+FP{i,j}(2));
                elseif ~isnan(FP{i,j}(1))
                    relFP(i,j) = -log(1+FP{i,j}(1));
                else
                    error('check');
                end
            end
        end
    end
    rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
    relFP = relFP(~rmInd,:);
    valid_rxns = targetRxns(~rmInd);
    % minus the min 
    relFP = relFP - min(relFP, [], 2);
    %Computing the correlation between a reaction expression and measured growth rate
    r=[];
    p_r=[];
    deltaminmax = [];
    testedRxn = {};
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat_normalized(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
            %Correcting for multiple hypothesis using FDR and significance level of
            deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
        end
    end
    
    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    
    n_corr_COMPASS(zz) = sum(fdr_r<0.05 & r >0);

end

end
