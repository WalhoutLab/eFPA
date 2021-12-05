%% About
% caculate the converted distance matrix that was used to convert weighted
% distance parameter to the real distance boudary of information integrated
% for every ROI target
%% 1. load the model and prepare the model
addpath('./../scripts/')
model = loadYeatModel();
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

proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);
%% flux table 
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;
% normalize flux unit
fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
%% make the converted dist mat
load('output/Titration_relativeExp_wtdDist_expDecay_FineGrained.mat');

cvtDistMat = nan(length(targetRxns),length(n2));
maxDist = max(distMat_normal(~isinf(distMat_normal)));
for zz = 1:length(n2)    
    FP = FP_collection_2{zz};
    %% convert the distance to real distance 
    for i = 1:length(targetRxns)
        if  mean(fluxMat(i,:)) > 0 
            if ~isnan(FP{i,end}(1))
                rxnID = [targetRxns{i},'_f'];
            elseif ~isnan(FP{i,end}(2))
                rxnID = [targetRxns{i},'_r'];
            else
                warning('check');
            end
        else
            if ~isnan(FP{i,end}(2))
                rxnID = [targetRxns{i},'_r'];
            elseif ~isnan(FP{i,end}(1))
                rxnID = [targetRxns{i},'_f'];
            else
                error('check');
            end
        end
        if any(strcmp(labels,rxnID))
            rxnSet = labels(distMat_wtd(strcmp(labels,rxnID),:) <= n2(zz));% these reactions were (expected to be) scaned through in the integration 
            cvtDistMat(i,zz) = max(distMat_normal(strcmp(labels,rxnID),ismember(labels,rxnSet))); % the maximum distance for those scaned rxns
        else
            cvtDistMat(i,zz) = maxDist;
        end
    end
end
save('output/cvtDistMat.mat','cvtDistMat');




