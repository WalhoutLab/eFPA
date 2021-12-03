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

%% flux table 
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;

fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
%% load merged levels
load('output/normalizedLevels_partialExcluded.mat');
%% first test one distance boundary 
n = 4;
intExp = nan(length(rxnLabel), length(conditions));
for i = 1:length(rxnLabel)
    targetRxn = rxnLabel{i};
    for j = 1:length(conditions)
        [level_f, level_r] = sumExp(targetRxn, n, labels, distMat_wtd, valid_rxns_pro_perPro, normalizedLevel_pro_perPro(:,j));
        if  mean(fluxMat(i,:)) > 0 
            if level_f ~= 0
                intExp(i,j) = level_f;
            else
                intExp(i,j) = level_r;
            end
        else
            if level_r ~= 0
                intExp(i,j) = level_r;
            else
                intExp(i,j) = level_f;
            end
        end
    end
end
%% evaluate it
rmInd = all(intExp == intExp(:,1),2);
intExp = intExp(~rmInd,:);
valid_rxns = rxnLabel(~rmInd);
%Computing the correlation between a reaction expression and measured growth rate
r=[];
p_r=[];
deltaminmax = [];
testedRxn = {};
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        [r(end+1),p_r(end+1)] = corr(intExp(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        deltaminmax(end+1) = max(intExp(strcmp(valid_rxns,rxnLabel{j}),:)) - min(intExp(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end

fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

%% titrate the boundaries and do the heatmap
% weighted distance
dorders = 0:0.5:40;
targetRxns = rxnLabel;
rMat = zeros(length(targetRxns),length(dorders));
FDRmat = ones(length(targetRxns),length(dorders));
pMat = ones(length(targetRxns),length(dorders));
CVmat = ones(length(targetRxns),length(dorders));
for nn = 1: length(dorders)
    %%
    n = dorders(nn);
    intExp = nan(length(rxnLabel), length(conditions));
    for i = 1:length(rxnLabel)
        targetRxn = rxnLabel{i};
        for j = 1:length(conditions)
            [level_f, level_r] = sumExp(targetRxn, n, labels, distMat_wtd, valid_rxns_pro_perPro, normalizedLevel_pro_perPro(:,j));
            if  mean(fluxMat(i,:)) > 0 
                if level_f ~= 0
                    intExp(i,j) = level_f;
                else
                    intExp(i,j) = level_r;
                end
            else
                if level_r ~= 0
                    intExp(i,j) = level_r;
                else
                    intExp(i,j) = level_f;
                end
            end
        end
    end

    
    rmInd = all(intExp == intExp(:,1),2);
    intExp = intExp(~rmInd,:);
    valid_rxns = rxnLabel(~rmInd);
    %Computing the correlation between a reaction expression and measured growth rate
    r=[];
    p_r=[];
    deltaminmax = [];
    testedRxn = {};
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat_normalized(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            [r(end+1),p_r(end+1)] = corr(intExp(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
            %Correcting for multiple hypothesis using FDR and significance level of
            deltaminmax(end+1) = max(intExp(strcmp(valid_rxns,rxnLabel{j}),:)) - min(intExp(strcmp(valid_rxns,rxnLabel{j}),:));
        end
    end

    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

    [A B] = ismember(testedRxn,targetRxns);
    rMat(B(A),nn) = r;
    FDRmat(B(A),nn) = fdr_r;
    CVmat(B(A),nn) = deltaminmax;
    pMat(B(A),nn) = p_r;
end
save('output/modeling_free_integration_weightedDist.mat', 'rMat','FDRmat','CVmat','pMat');
%% unweighted distance
% weighted distance
dorders = 0:0.5:40;
targetRxns = rxnLabel;
rMat = zeros(length(targetRxns),length(dorders));
FDRmat = ones(length(targetRxns),length(dorders));
pMat = ones(length(targetRxns),length(dorders));
CVmat = ones(length(targetRxns),length(dorders));
for nn = 1: length(dorders)
    %%
    n = dorders(nn);
    intExp = nan(length(rxnLabel), length(conditions));
    for i = 1:length(rxnLabel)
        targetRxn = rxnLabel{i};
        for j = 1:length(conditions)
            [level_f, level_r] = sumExp(targetRxn, n, labels, distMat_normal, valid_rxns_pro_perPro, normalizedLevel_pro_perPro(:,j));
            if  mean(fluxMat(i,:)) > 0 
                if level_f ~= 0
                    intExp(i,j) = level_f;
                else
                    intExp(i,j) = level_r;
                end
            else
                if level_r ~= 0
                    intExp(i,j) = level_r;
                else
                    intExp(i,j) = level_f;
                end
            end
        end
    end

    
    rmInd = all(intExp == intExp(:,1),2);
    intExp = intExp(~rmInd,:);
    valid_rxns = rxnLabel(~rmInd);
    %Computing the correlation between a reaction expression and measured growth rate
    r=[];
    p_r=[];
    deltaminmax = [];
    testedRxn = {};
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat_normalized(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            [r(end+1),p_r(end+1)] = corr(intExp(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
            %Correcting for multiple hypothesis using FDR and significance level of
            deltaminmax(end+1) = max(intExp(strcmp(valid_rxns,rxnLabel{j}),:)) - min(intExp(strcmp(valid_rxns,rxnLabel{j}),:));
        end
    end

    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

    [A B] = ismember(testedRxn,targetRxns);
    rMat(B(A),nn) = r;
    FDRmat(B(A),nn) = fdr_r;
    CVmat(B(A),nn) = deltaminmax;
    pMat(B(A),nn) = p_r;
    nn
end
save('output/modeling_free_integration_oriDist.mat', 'rMat','FDRmat','CVmat','pMat');
%% first plot the overall prediction 
dorders = 0:0.5:40;
load('output/modeling_free_integration_oriDist.mat')
addpath ./../6.humanTissue_MERGE/PlotPub/lib
baseline = 46;
N_sigCorr_highCV = [];
for nn = 1: length(dorders)
    N_sigCorr_highCV(1,nn) = sum(rMat((FDRmat(:,nn)<0.05 & CVmat(:,nn) > 0),nn)>0); % we dont put the 0.2 CV cutoff to make the control even more predictive
end
dorders_ori = dorders;
N_sigCorr_highCV_ori = N_sigCorr_highCV;

load('output/modeling_free_integration_weightedDist.mat')
N_sigCorr_highCV = [];
for nn = 1: length(dorders)
    N_sigCorr_highCV(1,nn) = sum(rMat((FDRmat(:,nn)<0.05 & CVmat(:,nn) > 0),nn)>0); % we dont put the 0.2 CV cutoff to make the control even more predictive
end
dorders_wtd = dorders;
N_sigCorr_highCV_wtd = N_sigCorr_highCV;

% load the control (FPA) 
load(['output/Titration_relativeExp_wtdDist_expDecay.mat'])
N_sigCorr = zeros(length(n2),length(base));
N_sigCorr_highCV = zeros(length(n2),length(base));
perc_sigCorr_in_highCV = zeros(length(n2),length(base));
for zz = 1:length(base)
    dorders = n2;
    rMat = zeros(length(targetRxns),length(dorders));
    for nn = 1: length(dorders)
        %%
        % nn = 10;
        %dorders(nn)
        FP = FP_collection_2{zz}{nn};
        relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
        for i = 1:size(FP,1)
            for j = 1:(size(FP,2)-1)
                if  mean(fluxMat(i,:)) > 0 
                    if ~isnan(FP{i,j}(1))
                        relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                    elseif ~isnan(FP{i,j}(2))
                        relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                    else
                        fprintf('base = %f, n = %d, i = %d, j = %d, is a failed FPA\n',...
                            base(zz), dorders(nn),i,j);
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
        fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

        [A B] = ismember(testedRxn,targetRxns);
        rMat(B(A),nn) = r;
        N_sigCorr(nn,zz) = sum(r(fdr_r<0.05)>0);
        N_sigCorr_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0);
        perc_sigCorr_in_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0) / sum(deltaminmax > 0.2);
    end
end

dorders_wtdDist_expDecay100 = dorders;% 1-dorders/max(dorders);
N_sigCorr_wtdDist_expDecay100 = N_sigCorr(:,7);% base =100
N_sigCorr_highCV_wtdDist_expDecay100 = N_sigCorr_highCV(:,7);
perc_sigCorr_in_highCV_wtdDist_expDecay100 = perc_sigCorr_in_highCV(:,7);

dorders_wtdDist_expDecay2 = dorders;% 1-dorders/max(dorders);
N_sigCorr_wtdDist_expDecay2 = N_sigCorr(:,2);% base =2
N_sigCorr_highCV_wtdDist_expDecay2 = N_sigCorr_highCV(:,2);
perc_sigCorr_in_highCV_wtdDist_expDecay2 = perc_sigCorr_in_highCV(:,2);

figure;
hold on
plot(dorders_wtdDist_expDecay2(1:end),N_sigCorr_highCV_wtdDist_expDecay2(1:end) ,'o-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
plot(dorders_wtdDist_expDecay100(1:end),N_sigCorr_highCV_wtdDist_expDecay100(1:end) ,'o-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
plot(dorders_wtd(1:end),N_sigCorr_highCV_wtd(1:end) ,'o-','LineWidth',2,'Color','#0072BD','MarkerSize', 3)
plot(dorders_ori(1:end),N_sigCorr_highCV_ori(1:end) ,'o-','LineWidth',2,'MarkerSize', 3)
% plot(dorders_expDecay100(1:13),N_sigCorr_highCV_expDecay100(1:13) ,'.--','LineWidth',2)
% plot(dorders_expDecay2(1:13),N_sigCorr_highCV_expDecay2(1:13) ,'.--','LineWidth',2)
yline(baseline,'--','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Distance order/boundary');
ylabel('Number of significantly correlated reactions ');
ylim([20 80])
legend({'improved FPA (base = 2)','improved FPA (base = 100)','local expression average (weighted dist)','local expression average (naive dist)','target expression only'},'FontSize',7);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.LegendLoc = 'NorthEast';
plt.export(['figures/modeling_free_integration.pdf']);
%% plot the titration heatmap
load('output/modeling_free_integration.mat')
dorders = 0:0.5:40;

pMat_valid = [pMat];
rMat_valid = [rMat];
CVmat_valid = [CVmat];
% FDRmat_valid(CVmat<0.1) = 1;
% rMat_valid(CVmat<0.1) = 0;

keep = any(pMat_valid < 0.05 & rMat_valid > 0 & CVmat_valid > 0.2,2);

rMat_valid = rMat_valid(keep,:);
pMat_valid = pMat_valid(keep,:);
targetRxns_valid = targetRxns(keep);
CVmat_valid = CVmat_valid(keep,:);


% rMat_valid_normalized = -log10(FDRmat_valid) .* sign(rMat_valid);
rMat_valid_normalized = rMat_valid ./ max(rMat_valid,[],2); % relative r 
%rMat_valid ./ max(rMat_valid,[],2);% normalize(rMat_valid,2,'range');
% rMat_valid_normalized = rMat_valid_normalized ./ max(rMat_valid_normalized,[],2);



distMethod = 'euclidean';
IDs = [strsplit(num2str(dorders))];

cgo=clustergram(rMat_valid_normalized(:,1:end),'RowLabels',targetRxns_valid,'ColumnLabels',IDs(1:end),'RowPDist',distMethod,'Cluster', 'Column');
c=get(cgo,'ColorMap');
n = 100;
tmp = [zeros(n,1), linspace(1,0,n)',zeros(n,1)];
tmp = tmp(2:end,:);
cpr=[tmp; linspace(0,1,n)',zeros(n,1),zeros(n,1)
    ];
% cpr = cpr(199:-1:1,:);
set(cgo,'ColorMap',cpr);
set(cgo,'Symmetric',false);

set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 12)

