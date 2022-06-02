%% About
% FPA evaluation. systematically compare the FPA prediction with no
% integration and evalute the distance boundaries
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/
addpath ./../../../PlotPub/lib/
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
distMat = distMat_min;
% load the special penalties 
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
%% load expression levels
load('output/normalizedLevels_partialExcluded.mat');
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
%% the second evaluation set (156 rxns with expression measured)
mySet = intersect(valid_rxns_pro_perPro, rxnLabel);
%% gather the finished randomization files
allFiles = dir('output/randomization_simpleDecay_seed_controlled/seed_*/rand_*_default_FPA.mat');
names = {allFiles.name};
folders = {allFiles.folder};
folders = folders(~strcmp(names,'rand_1_default_FPA.mat'));
names = names(~strcmp(names,'rand_1_default_FPA.mat'));
for i = 1:length(names)
    path_defaultFPA{i} = [folders{i},'/',names{i}];
end

allFiles = dir('output/randomization_simpleDecay_seed_controlled/seed_*/rand_*_flexi_FPA.mat');
names = {allFiles.name};
folders = {allFiles.folder};
folders = folders(~strcmp(names,'rand_1_flexi_FPA.mat'));
names = names(~strcmp(names,'rand_1_flexi_FPA.mat'));
for i = 1:length(names)
    path_flexiFPA{i} = [folders{i},'/',names{i}];
end
% only read completed data
rmInd = [];
for i = 1:length(path_defaultFPA)
    mypath = path_defaultFPA{i};
    mypath = regexprep(mypath, 'default_FPA.mat$','flexi_FPA.mat');
    if ~any(strcmp(mypath, path_flexiFPA))
        rmInd = [rmInd; i];
    end
end
path_defaultFPA(rmInd) = [];
% we only keep the first 1000 randomization
path_defaultFPA = path_defaultFPA(1:1000);
path_flexiFPA = path_flexiFPA(1:1000);
%% evaluate FPA randomization for default FPA
n2 = 0:0.5:40;
N_corr_232 = [];
N_corr_156 = [];
for zz = 1:length(path_defaultFPA)
    %% correlation for 1-d titration
    load(path_defaultFPA{zz}) 
    FP = FP_collection_2{1}{1};
    relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
    for i = 1:size(FP,1)
        for j = 1:(size(FP,2)-1)
            if  mean(fluxMat(i,:)) > 0 
                if ~isnan(FP{i,j}(1))
                    relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                elseif ~isnan(FP{i,j}(2))
                    relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                else
                    warning('check');
                end
            else
                if ~isnan(FP{i,j}(2))
                    relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                elseif ~isnan(FP{i,j}(1))
                    relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                else
                    error('check');
                end
            end
        end
    end
    rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
    relFP = relFP(~rmInd,:);
    valid_rxns = targetRxns(~rmInd);
    % Computing the correlation
    r=[];
    p_r=[];
    deltaminmax = [];
    testedRxn = {};
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat_normalized(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
            deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
        end
    end

    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0));


    N_corr_232(zz) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0);
    N_corr_156(zz) = sum(fdr_r<0.05 & deltaminmax > 0.2 & r > 0 & ismember(testedRxn, mySet));
end

%% plot the benchmark by 232 (all) reactions
baseline1 = 73; % number of correlated rxns by ROI expression only
baseline2 = 52; % number of correlated rxns by ROI expression only

figure;
hold on
histogram(N_corr_232)
xline(baseline1,'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Number of significantly correlated reactions');
ylabel('Number of permutation');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_defaultFPA_232.pdf']);

figure;
hold on
histogram(N_corr_156)
xline(baseline2,'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Number of significantly correlated reactions');
ylabel('Number of permutation');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_defaultFPA_156.pdf']);
%% analyze the randomization of titration
n2 = 0:0.5:40;
N_corr_232 = [];
N_corr_156 = [];
N_corr_232_cum = []; % cummulative predicted rxns as distance bound goes up
N_corr_156_cum = [];
Nmat_corr_232 = [];
Nmat_corr_156 = [];
rMax = [];
for zz = 1:length(path_flexiFPA)
    load(path_flexiFPA{zz}) 
    dorders = n2;
    rMat = zeros(length(targetRxns),length(dorders));
    FDRmat = ones(length(targetRxns),length(dorders));
    pMat = ones(length(targetRxns),length(dorders));
    CVmat = zeros(length(targetRxns),length(dorders));
    FP_collection_2 = FP_collection_2{1};
    for nn = 1: length(dorders)
        %%
        FP = FP_collection_2{nn};
        relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
        for i = 1:size(FP,1)
            for j = 1:(size(FP,2)-1)
                if  mean(fluxMat(i,:)) > 0 
                    if ~isnan(FP{i,j}(1))
                        relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                    elseif ~isnan(FP{i,j}(2))
                        relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                    else
                        warning('check');
                    end
                else
                    if ~isnan(FP{i,j}(2))
                        relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                    elseif ~isnan(FP{i,j}(1))
                        relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                    else
                        error('check');
                    end
                end
            end
        end
        rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
        relFP = relFP(~rmInd,:);
        valid_rxns = targetRxns(~rmInd);
        %Computing the correlation 
        r=[];
        p_r=[];
        deltaminmax = [];
        testedRxn = {};
        for j = 1:length(rxnLabel)
            fluxMeasure = fluxMat_normalized(j,:);
            if any(strcmp(valid_rxns,rxnLabel{j}))
                testedRxn(end+1) = rxnLabel(j);
                [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
                deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
            end
        end

        fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
        fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0));

        [A B] = ismember(testedRxn,targetRxns);
        rMat(B(A),nn) = r;
        FDRmat(B(A),nn) = fdr_r;
        CVmat(B(A),nn) = deltaminmax;
        pMat(B(A),nn) = p_r;
        
    end
    N_corr_232(zz) = sum(any(FDRmat < 0.05 & rMat > 0 & CVmat > 0.2,2));
    N_corr_156(zz) = sum(ismember(targetRxns, mySet) & any(FDRmat < 0.05 & rMat > 0 & CVmat > 0.2,2));
    Nmat_corr_232(zz,:) = sum(FDRmat < 0.05 & rMat > 0 & CVmat > 0.2,1);
    Nmat_corr_156(zz,:) = sum(FDRmat(ismember(targetRxns, mySet),:) < 0.05 & rMat(ismember(targetRxns, mySet),:) > 0 & CVmat(ismember(targetRxns, mySet),:) > 0.2,1);
    for k = 1:length(n2)
        N_corr_232_cum(zz,k) = sum(any(FDRmat(:,1:k) < 0.05 & rMat(:,1:k) > 0 & CVmat(:,1:k) > 0.2,2));
        N_corr_156_cum(zz,k) = sum(ismember(targetRxns, mySet) & any(FDRmat(:,1:k) < 0.05 & rMat(:,1:k) > 0 & CVmat(:,1:k) > 0.2,2));
    end
    rMax(zz,:) = max(rMat,[],2);
end

%% the control
load('output/randomization_simpleDecay/seed_1126/rand_1_flexi_FPA.mat') 
dorders = n2;
rMat = zeros(length(targetRxns),length(dorders));
FDRmat = ones(length(targetRxns),length(dorders));
pMat = ones(length(targetRxns),length(dorders));
CVmat = zeros(length(targetRxns),length(dorders));
FP_collection_2 = FP_collection_2{1};
for nn = 1: length(dorders)
    %%
    FP = FP_collection_2{nn};
    relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
    for i = 1:size(FP,1)
        for j = 1:(size(FP,2)-1)
            if  mean(fluxMat(i,:)) > 0 
                if ~isnan(FP{i,j}(1))
                    relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                elseif ~isnan(FP{i,j}(2))
                    relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                else
                    warning('check');
                end
            else
                if ~isnan(FP{i,j}(2))
                    relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                elseif ~isnan(FP{i,j}(1))
                    relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                else
                    error('check');
                end
            end
        end
    end
    rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
    relFP = relFP(~rmInd,:);
    valid_rxns = targetRxns(~rmInd);
    %Computing the correlation 
    r=[];
    p_r=[];
    deltaminmax = [];
    testedRxn = {};
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat_normalized(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
            deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
        end
    end

    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0));

    [A B] = ismember(testedRxn,targetRxns);
    rMat(B(A),nn) = r;
    FDRmat(B(A),nn) = fdr_r;
    CVmat(B(A),nn) = deltaminmax;
    pMat(B(A),nn) = p_r;

end
N_corr_232_obs = sum(any(FDRmat < 0.05 & rMat > 0 & CVmat > 0.2,2));
N_corr_156_obs = sum(ismember(targetRxns, mySet) & any(FDRmat < 0.05 & rMat > 0 & CVmat > 0.2,2));
Nmat_corr_232_obs = sum(FDRmat < 0.05 & rMat > 0 & CVmat > 0.2,1);
Nmat_corr_156_obs = sum(FDRmat(ismember(targetRxns, mySet),:) < 0.05 & rMat(ismember(targetRxns, mySet),:) > 0 & CVmat(ismember(targetRxns, mySet),:) > 0.2,1);
for k = 1:length(n2)
    N_corr_232_cum_obs(k) = sum(any(FDRmat(:,1:k) < 0.05 & rMat(:,1:k) > 0 & CVmat(:,1:k) > 0.2,2));
    N_corr_156_cum_obs(k) = sum(ismember(targetRxns, mySet) & any(FDRmat(:,1:k) < 0.05 & rMat(:,1:k) > 0 & CVmat(:,1:k) > 0.2,2));
end
rMax_obs = max(rMat,[],2);
%% plot
figure;
hold on
histogram(N_corr_232, 'NumBins',30)
xline(N_corr_232_obs,'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Number of significantly correlated reactions');
ylabel('Number of permutation');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.XLim = [15, 120]
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_232.pdf']);

figure;
hold on
histogram(N_corr_156)
xline(N_corr_156_obs,'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Number of significantly correlated reactions');
ylabel('Number of permutation');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_156.pdf']);


figure;
hold on
%plot(n2,Nmat_corr_156_obs ,'-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
for i = 1:length(path_flexiFPA)
    plot(n2,Nmat_corr_232(i,:) ,'-','LineWidth',2,'MarkerSize', 3,'Color',[0,0,0,0.2])
    %plot(n2,Nmat_corr_156(i,:) ,'-k','LineWidth',2,'MarkerSize', 3)
end
plot(n2,Nmat_corr_232_obs,'-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
xlabel('Distance boundary');
ylabel('Number of significantly correlated reactions ');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85*2, 2.35*2];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_titration.pdf']);

figure;
hold on
%plot(n2,N_corr_156_cum_obs ,'-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
for i = 1:length(path_flexiFPA)
    plot(n2,N_corr_232_cum(i,:) ,'-','LineWidth',2,'MarkerSize', 3,'Color',[0,0,0,0.2])
    %plot(n2,N_corr_156_cum(i,:) ,'-k','LineWidth',2,'MarkerSize', 3)
end
plot(n2,N_corr_232_cum_obs,'-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
xlabel('Distance boundary');
ylabel('Cummulative number of significantly correlated reactions ');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85*2, 2.35*2];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_titration_cummulative.pdf']);

% note: in observed data, 94/101 reactions were predicted with distance
% boundary less or equal to 15, which is biologically meaningful long
% pathway, but randomised data predicted few until distance boundary is
% greater than 15.

% the two examples
% r_0362 - glycan
figure;
hold on
histogram(rMax(:,strcmp(targetRxns,'r_0362')), 'NumBins',30)
xline(rMax_obs(strcmp(targetRxns,'r_0362')),'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Max PCC in boundary titration');
ylabel('Number of permutation');
title(['p < ',num2str((sum(rMax(:,strcmp(targetRxns,'r_0362')) >= rMax_obs(strcmp(targetRxns,'r_0362')))+1) ./ (size(rMax,1)+1),2)])
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_glycan_example_r_0362.pdf']);

% r_0485 - gst
figure;
hold on
histogram(rMax(:,strcmp(targetRxns,'r_0485')), 'NumBins',30)
xline(rMax_obs(strcmp(targetRxns,'r_0485')),'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Max PCC in boundary titration');
ylabel('Number of permutation');
title(['p < ',num2str((sum(rMax(:,strcmp(targetRxns,'r_0485')) >= rMax_obs(strcmp(targetRxns,'r_0485')))+1) ./ (size(rMax,1)+1),2)])
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_glycan_example_r_0485.pdf']);
%% the increase compared with initial status
figure;
hold on
%plot(n2,Nmat_corr_156_obs ,'-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
for i = 1:length(path_flexiFPA)
    plot(n2,(Nmat_corr_232(i,:)-Nmat_corr_232(i,1)) ./ Nmat_corr_232(i,1) ,'-','LineWidth',2,'MarkerSize', 3,'Color',[0,0,0,0.2])
    %plot(n2,Nmat_corr_156(i,:) ,'-k','LineWidth',2,'MarkerSize', 3)
end
plot(n2,(Nmat_corr_232_obs - Nmat_corr_232_obs(1)) ./ Nmat_corr_232_obs(1),'-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
xlabel('Distance boundary');
ylabel('Proportional increase (% of the starting point)');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85*2, 2.35*2];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_proportional_benifit.pdf']);

figure;
hold on
%plot(n2,N_corr_156_cum_obs ,'-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
for i = 1:length(path_flexiFPA)
    plot(n2,(N_corr_232_cum(i,:) - N_corr_232_cum(i,1))./ N_corr_232_cum(i,1) ,'-','LineWidth',2,'MarkerSize', 3,'Color',[0,0,0,0.2])
    %plot(n2,N_corr_156_cum(i,:) ,'-k','LineWidth',2,'MarkerSize', 3)
end
plot(n2,(N_corr_232_cum_obs - N_corr_232_cum_obs(1)) ./ N_corr_232_cum_obs(1),'-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
xlabel('Distance boundary');
ylabel('Cummulative proportional increase (% of the starting point)');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85*2, 2.35*2];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_cummulative_proportional_benifit.pdf']);


figure;
hold on
%plot(n2,Nmat_corr_156_obs ,'-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
for i = 1:length(path_flexiFPA)
    plot(n2,(Nmat_corr_232(i,:)-Nmat_corr_232(i,1)) ,'-','LineWidth',2,'MarkerSize', 3,'Color',[0,0,0,0.2])
    %plot(n2,Nmat_corr_156(i,:) ,'-k','LineWidth',2,'MarkerSize', 3)
end
plot(n2,(Nmat_corr_232_obs - Nmat_corr_232_obs(1)) ,'-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
xlabel('Distance boundary');
ylabel('Absolute increase (# of reactions cmpd w/ the starting point)');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85*2, 2.35*2];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_absolute_benifit.pdf']);

figure;
hold on
%plot(n2,N_corr_156_cum_obs ,'-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
for i = 1:length(path_flexiFPA)
    plot(n2,(N_corr_232_cum(i,:) - N_corr_232_cum(i,1)),'-','LineWidth',2,'MarkerSize', 3,'Color',[0,0,0,0.2])
    %plot(n2,N_corr_156_cum(i,:) ,'-k','LineWidth',2,'MarkerSize', 3)
end
plot(n2,(N_corr_232_cum_obs - N_corr_232_cum_obs(1)),'-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
xlabel('Distance boundary');
ylabel('Cummulative absolute increase (# of reactions cmpd w/ the starting point)');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85*2, 2.35*2];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_flexiFPA_cummulative_absolute_benifit.pdf']);


figure;
hold on
histogram(N_corr_232)
xline(baseline1,'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Number of significantly correlated reactions');
ylabel('Number of permutation');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_defaultFPA_232.pdf']);
%% get the empirical p-value of each reaction's prediction (reconciled or not by optimal bound FPA)
% we only consider the PCC for simplicity for now. maybe we will use the
% same criteria as we define the reconcilitation (PCC > 0, FDR < 0.05 range
% > 0.2) that disable the histogram analysis but provided real p-value
for i = 1:size(rMax,2)
    pVals(i) = (sum(rMax(:,i) >= rMax_obs(i))+1) ./ (size(rMax,1)+1);
end
sum(pVals < 0.05)
sum(Vals < 0.1)
sum(pVals < 0.25)
% we stop overinterpreting this...


%% optional
% first manually rerun the section "evaluate FPA randomization for default FPA" to reset
% variables

% also plot the net increase
figure;
hold on
histogram(N_corr_232 - N_corr_232_cum(:,1)')
xline(baseline1 - 46,'-','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Number of newly predicted reactions');
ylabel('Number of permutation');
title(['p < ',num2str((sum( N_corr_232 - N_corr_232_cum(:,1)' >= baseline1 - 46 )+1) ./ (length(N_corr_232)+1),2)])
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.export(['figures/NoTrack_randomization_defaultFPA_232_absolute_benefit.pdf']);

