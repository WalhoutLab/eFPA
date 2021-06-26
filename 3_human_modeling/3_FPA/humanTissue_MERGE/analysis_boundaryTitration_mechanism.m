%% load data
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;
%% caculate the mean maxIRS for each reaction 
heatTbl = readtable('output/heatmapTbl.csv');
targetRxns = heatTbl.Row;
load('output/Titration_relativeExp_wtdDist_expDecay_FineGrained.mat');

mmIRS = nan(length(targetRxns),length(n2));
for zz = 1:length(n2)
    %% load the FPA table set 
    load(['fine_titration/base_1000_n_', num2str(n2(zz)),'.mat']);
    
    %% clean up the predictions (only keep the relavant directions)
    FPsol = cell(size(FluxPotential_solutions,1),size(FluxPotential_solutions,2));%we choose one from f and r as a prediction
    for i = 1:size(FluxPotential_solutions,1)
        for j = 1:(size(FluxPotential_solutions,2))
            if  mean(fluxMat(i,:)) > 0 
                if ~isempty(FluxPotential_solutions{i,j}(1))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(1);
                elseif ~isempty(FluxPotential_solutions{i,j}(2))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(2);
                else
                    fprintf('i = %d, j = %d, is a failed FPA\n',...
                        i,j);
                    warning('check');
                end
            else
                if ~isempty(FluxPotential_solutions{i,j}(2))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(2);
                elseif ~isempty(FluxPotential_solutions{i,j}(1))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(1);
                else
                    error('check');
                end
            end
        end
    end
    %% check distribution of influential boudary
    % we define the Influential Reaction Set (IRS) as the top flux allowance
    % contributors that have a accumulative flux allowance contribution of 95%
    % for each condition, we consider three metrics of the natural distance of
    % all IRS reactions: mean, median, maximum 
    meanIRS = nan(size(FPsol,1),size(FPsol,2));
    medianIRS = nan(size(FPsol,1),size(FPsol,2));
    maximumIRS = nan(size(FPsol,1),size(FPsol,2));

    for i = 1:size(FPsol,1)
        for j = 1:size(FPsol,2)
            res = FPsol{i,j};
            mytbl = struct();
            mytbl.labels = res.labels;
            mytbl.fluxAllowance = res.full .* res.weight;
            mytbl.UnwDist = res.UnwDist';
            mytbl = struct2table(mytbl);
            mytbl = sortrows(mytbl,2,'descend');
            IRS = mytbl.UnwDist(1);
            accFA = mytbl.fluxAllowance(1);
            for k = 2:size(mytbl,1)
               accFA = accFA + mytbl.fluxAllowance(k);
               if accFA < 0.95
                   IRS = [IRS,mytbl.UnwDist(k)];
               else
                   break;
               end
            end
            meanIRS(i,j) = mean(IRS);
            medianIRS(i,j) = median(IRS);
            maximumIRS(i,j) = max(IRS);
        end
    end

    % ==> maxIRS is the most stable metric, and it is also reasonable
    % so we use the mean of maxIRS as the distance dictator of the mechanism of
    % each reaction. 
    mmIRS(:,zz) = mean(maximumIRS,2);
    zz
end
save('output/mmIRS.mat','mmIRS');

%% same for base2 n6
%% caculate the mean maxIRS for each reaction 
load(['base2_n6_eng/base_2_n_6.mat']);
% clean up the predictions (only keep the relavant directions)
FPsol = cell(size(FluxPotential_solutions,1),size(FluxPotential_solutions,2));%we choose one from f and r as a prediction
for i = 1:size(FluxPotential_solutions,1)
    for j = 1:(size(FluxPotential_solutions,2))
        if  mean(fluxMat(i,:)) > 0 
            if ~isempty(FluxPotential_solutions{i,j}(1))
                FPsol(i,j) = FluxPotential_solutions{i,j}(1);
            elseif ~isempty(FluxPotential_solutions{i,j}(2))
                FPsol(i,j) = FluxPotential_solutions{i,j}(2);
            else
                fprintf('i = %d, j = %d, is a failed FPA\n',...
                    i,j);
                warning('check');
            end
        else
            if ~isempty(FluxPotential_solutions{i,j}(2))
                FPsol(i,j) = FluxPotential_solutions{i,j}(2);
            elseif ~isempty(FluxPotential_solutions{i,j}(1))
                FPsol(i,j) = FluxPotential_solutions{i,j}(1);
            else
                error('check');
            end
        end
    end
end
% check distribution of influential boudary
% we define the Influential Reaction Set (IRS) as the top flux allowance
% contributors that have a accumulative flux allowance contribution of 95%
% for each condition, we consider three metrics of the natural distance of
% all IRS reactions: mean, median, maximum 
meanIRS = nan(size(FPsol,1),size(FPsol,2));
medianIRS = nan(size(FPsol,1),size(FPsol,2));
maximumIRS = nan(size(FPsol,1),size(FPsol,2));

for i = 1:size(FPsol,1)
    for j = 1:size(FPsol,2)
        res = FPsol{i,j};
        mytbl = struct();
        mytbl.labels = res.labels;
        mytbl.fluxAllowance = res.full .* res.weight;
        mytbl.UnwDist = res.UnwDist';
        mytbl = struct2table(mytbl);
        mytbl = sortrows(mytbl,2,'descend');
        IRS = mytbl.UnwDist(1);
        accFA = mytbl.fluxAllowance(1);
        for k = 2:size(mytbl,1)
           accFA = accFA + mytbl.fluxAllowance(k);
           if accFA < 0.95
               IRS = [IRS,mytbl.UnwDist(k)];
           else
               break;
           end
        end
        meanIRS(i,j) = mean(IRS);
        medianIRS(i,j) = median(IRS);
        maximumIRS(i,j) = max(IRS);
    end
end

mmIRS_b2n6 = mean(maximumIRS,2);

save('output/mmIRS_base2n6.mat','mmIRS_b2n6');

%% manual inspection of a table
% res = FluxPotential_solutions{10,1}{1};
% mytbl = struct();
% mytbl.labels = res.labels;
% mytbl.full = res.full;
% mytbl.fluxAllowance = mytbl.full .* res.weight;
% mytbl.weight = res.weight;
% mytbl.penalty = res.penalty;
% mytbl.pDist = res.pDist';
% mytbl.Dist = res.Dist';
% mytbl.UnwDist = res.UnwDist';
% mytbl = struct2table(mytbl);
% mytbl = sortrows(mytbl,3,'descend');

%% caculate the converted distance matrix 
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
            rxnSet = labels(distMat_wtd(strcmp(labels,rxnID),:) <= n2(zz));
            cvtDistMat(i,zz) = max(distMat_normal(strcmp(labels,rxnID),ismember(labels,rxnSet)));
        else
            cvtDistMat(i,zz) = maxDist;
        end
    end
end
save('output/cvtDistMat.mat','cvtDistMat');




