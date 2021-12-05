%% About
% preparing the tables for the distance titration heatmap of PCC. We also
% made the conversion of weighted metabolic distance of the distance
% boudary parameter to interpretable real distance in the plot 
%%
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/
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
%% load the PCC by ROI expression only
normalizedLevel = normalizedLevel_pro_perPro;
valid_rxns = valid_rxns_pro_perPro ;
%Computing the correlation
rho=[];
p_rho=[];
r=[];
p_r=[];
testedRxn = {};
rxnLabel = fluxTbl.Model_Reaction_ID;
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        [r(end+1),p_r(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
    end
end
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
control_r = r;
control_fdr = fdr_r;
control_testedRxns = testedRxn;
%% load the default improved FPA (local integration) - base 2 boundary 6
load(['output/Titration_relativeExp_wtdDist_expDecay.mat'])
dorders = n2;
nn = 9;

FP = FP_collection_2{2}{nn};
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
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

control_r_2 = r;
control_fdr_2 = fdr_r;
control_testedRxns_2 = testedRxn;
control_deltaminmax_2 = deltaminmax;
%% load the distance bound titration and make the heatmap 
load(['output/Titration_relativeExp_wtdDist_expDecay_FineGrained.mat'])
%% calculate correlation
dorders = n2;
rMat = zeros(length(targetRxns),length(dorders));
FDRmat = ones(length(targetRxns),length(dorders));
pMat = ones(length(targetRxns),length(dorders));
CVmat = ones(length(targetRxns),length(dorders));
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

%% merge the results and write out tables
% ROI expression only
[A B] = ismember(targetRxns,control_testedRxns);
ctrVp = ones(length(targetRxns),1);
ctrVr = zeros(length(targetRxns),1);
ctrVdeltaminmax = zeros(length(targetRxns),1);
ctrVp(A) = control_fdr(B(A)); % actually fdr
ctrVr(A) = control_r(B(A));
ctrVdeltaminmax(A) = 1; % we dont filter the deltaminmax of expression only result

% default FPA
[A B] = ismember(targetRxns,control_testedRxns_2);
ctrVp2 = ones(length(targetRxns),1);
ctrVr2 = zeros(length(targetRxns),1);
ctrVdeltaminmax2 = zeros(length(targetRxns),1);
ctrVp2(A) = control_fdr_2(B(A));
ctrVr2(A) = control_r_2(B(A));
ctrVdeltaminmax2(A) = control_deltaminmax_2(B(A));

% merged
pMat_valid = [ctrVp2, ctrVp, FDRmat];
rMat_valid = [ctrVr2, ctrVr, rMat];
deltaminmaxmat_valid = [ctrVdeltaminmax2, ctrVdeltaminmax, CVmat];
t = array2table(rMat_valid);
t.Properties.RowNames = targetRxns;
writetable(t,'output/PCC_titration_all.csv','WriteRowNames',1);
%% plot the histogram of max PCC (optimal-boundary integration PCC)
r_max = [];
r_max_sig = [];
for i = 1:size(rMat,1)
    if all(rMat(i,:) == 0)
        r_max(i) = NaN;
        r_max_sig(i) = NaN;
    else
        r_pass = rMat(i,CVmat(i,:) > 0.2); % we only consider valid predictions (range > 0.2)
        fdr_pass = FDRmat(i,CVmat(i,:) > 0.2);
        if(length(r_pass)> 0)
            r_max(i) = max(r_pass);
            r_max_sig(i) = min(fdr_pass(r_pass == max(r_pass)));
        else
            r_max(i) = max(rMat(i,:));
            r_max_sig(i) = 1;% this prediction will have a FDR of 1 so that it will never be considered as significant predcition
        end
    end
end
figure;
hold on
histogram(r_max,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(r_max(r_max_sig<0.05 & r_max >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
hold off
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export('figures/optimal_boundary_FPA_flux_correlation_pearson.pdf');
%% filter the PCC matrix and only keep the predicted reactions 
keep = any(pMat_valid < 0.05 & rMat_valid > 0 & deltaminmaxmat_valid > 0.2,2);

rMat_valid = rMat_valid(keep,:);
pMat_valid = pMat_valid(keep,:);
targetRxns_valid = targetRxns(keep);
deltaminmaxmat_valid = deltaminmaxmat_valid(keep,:);

% row-wise normalization
rMat_valid_normalized = rMat_valid ./ max(rMat_valid,[],2); % relative r 
rMat_valid_normalized_ori = rMat_valid_normalized;

% label the significantly correlated predictions 
sigMat = pMat_valid < 0.05 & rMat_valid > 0 & deltaminmaxmat_valid > 0.2;
%% write out the PCC matrix whose distance label was converted to real distance 
% bin the PCC according to converted distance 
dorders = 0:40;
rMat_valid_normalized = rMat_valid_normalized_ori;
load('output/cvtDistMat.mat');
cvtDistMat = cvtDistMat(keep,:);
rMat_valid_normalized_binned = nan(size(rMat_valid_normalized,1),length(dorders)+2);
sigMat_binned = nan(size(rMat_valid_normalized,1),length(dorders)+2);
for i = 1:(length(dorders)-1)
    lb = dorders(i);
    ub = dorders(i+1);
    for j = 1:size(rMat_valid_normalized_binned,1)
        pass = find(cvtDistMat(j,:) >= lb & cvtDistMat(j,:) < ub)+2;
        if any(pass)
            rMat_valid_normalized_binned(j,i+2) = max(rMat_valid_normalized(j,pass));
            sigMat_binned(j,i+2) = any(pMat_valid(j,pass) < 0.05 & rMat_valid(j,pass) > 0 & deltaminmaxmat_valid(j,pass) > 0.2);
        end
    end
end
rMat_valid_normalized_binned(:,1:2) = rMat_valid_normalized(:,1:2);
sigMat_binned(:,1:2) = sigMat(:,1:2);
% save data 
t = array2table(rMat_valid_normalized_binned);
t.Properties.RowNames = targetRxns_valid;
writetable(t,'output/relCorr_heatmapTbl_realDist.csv','WriteRowNames',1);

IDs = [{'base 2 - boundary 6','expression only'},strsplit(num2str(dorders))];
boundaries =  cell2table(IDs');
writetable(boundaries,'output/heatmapTbl_boundaries_realDist.csv');

t = array2table(sigMat_binned);
t.Properties.RowNames = targetRxns_valid;
writetable(t,'output/heatmapTbl_sigLabel_realDist.csv','WriteRowNames',1);