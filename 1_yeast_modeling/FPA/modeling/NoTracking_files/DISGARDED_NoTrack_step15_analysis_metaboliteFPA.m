%% About
% preparing the tables for the distance titration heatmap of PCC. We also
% made the conversion of weighted metabolic distance of the distance
% boudary parameter to interpretable real distance in the plot 
%%
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
%% metabolite table (relative metabolite in log2 scale)
% metTbl = readtable('./../input/YeastJoshua/originalDataTbl/metTbl.xlsx');
% % unify id 
% metTbl.Metabolite(strcmp(metTbl.Metabolite,'3-phosphoglycerate')) = {'3-phosphonato-D-glycerate(3-)'};
% metTbl.Metabolite(strcmp(metTbl.Metabolite,'alpha-D-ribose 1-phosphate')) = {'alpha-D-ribose 1-phosphate(2-)'};
% % use the met name in following analyses. 
% metTbl.Model_Metabolite_ID = strcat(metTbl.Metabolite,' [cytoplasm]');
% writetable(metTbl, './../input/YeastJoshua/originalDataTbl/metTbl.csv');
% id check
% we map metabolites by their names
% setdiff(metTbl.Metabolite,regexprep(model.metNames,' \[.+\]$',''))
% 3-phosphoglycerate ==> 3-phosphonato-D-glycerate(3-)
% alpha-D-ribose 1-phosphate ==> alpha-D-ribose 1-phosphate(2-)
metTbl = readtable('./../input/YeastJoshua/originalDataTbl/metTbl.csv');
% make the matched matrix
for i = 1: length(conditions)
    metMat(:,i) = metTbl.(conditions{i});
end
metLabel = metTbl.Model_Metabolite_ID;
% change to regular scale
metMat = 2.^metMat;
%% load the flux matrix
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
%% load the default improved FPA (local integration) - base 2 boundary 6
load('output/metaboliteFPA_default_FPA.mat')

FP = FP_collection_2{1}{1};
relFP_prod = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
relFP_cons = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
for i = 1:size(FP,1)
    for j = 1:(size(FP,2)-1)
        relFP_prod(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_cons(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
rmInd1 = all(isnan(relFP_prod),2) | all(relFP_prod == relFP_prod(:,1),2);
rmInd2 = all(isnan(relFP_cons),2) | all(relFP_cons == relFP_cons(:,1),2);
relFP_prod = relFP_prod(~rmInd1,:);
valid_mets_prod = regexprep(targetRxns(~rmInd1),'^NewMet_','');
relFP_cons = relFP_cons(~rmInd2,:);
valid_mets_cons= regexprep(targetRxns(~rmInd2),'^NewMet_','');

%Computing the correlation - production
r=[];
p_r=[];
deltaminmax = [];
testedMet = {};
for j = 1:length(metLabel)
    metMeasure = metMat(j,:);
    if any(strcmp(valid_mets_prod,metLabel{j}))
        testedMet(end+1) = metLabel(j);
        [r(end+1),p_r(end+1)] = corr(relFP_prod(strcmp(valid_mets_prod,metLabel{j}),:)',metMeasure','type','Pearson');
        deltaminmax(end+1) = max(relFP_prod(strcmp(valid_mets_prod,metLabel{j}),:)) - min(relFP_prod(strcmp(valid_mets_prod,metLabel{j}),:));
    end
end
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d mets give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));
localFPA_r_prod = r;
localFPA_fdr_prod = fdr_r;
localFPA_testedMets_prod = testedMet;
localFPA_deltaminmax_prod = deltaminmax;
% quick plot
% myset = testedMet(r > 0 & fdr_r < 0.05)';
% for j = 1:length(metLabel)
%     metMeasure = metMat(j,:);
%     if any(strcmp(myset,metLabel{j}))
%         lm = fitlm(relFP_prod(strcmp(valid_mets_prod,metLabel{j}),:)',metMeasure');
%         figure;
%         plot(lm)
%         title(metLabel{j})
%     end
% end

%Computing the correlation - consumption
r=[];
p_r=[];
deltaminmax = [];
testedMet = {};
for j = 1:length(metLabel)
    metMeasure = metMat(j,:);
    if any(strcmp(valid_mets_cons,metLabel{j}))
        testedMet(end+1) = metLabel(j);
        [r(end+1),p_r(end+1)] = corr(relFP_cons(strcmp(valid_mets_cons,metLabel{j}),:)',metMeasure','type','Pearson');
        deltaminmax(end+1) = max(relFP_cons(strcmp(valid_mets_cons,metLabel{j}),:)) - min(relFP_cons(strcmp(valid_mets_cons,metLabel{j}),:));
    end
end
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d mets give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));
localFPA_r_cons = r;
localFPA_fdr_cons = fdr_r;
localFPA_testedMets_cons = testedMet;
localFPA_deltaminmax_cons = deltaminmax;
allPredicted_localFPA = union(localFPA_testedMets_cons(localFPA_r_cons > 0 & localFPA_fdr_cons < 0.05)', ...
    localFPA_testedMets_prod(localFPA_r_prod > 0 & localFPA_fdr_prod < 0.05)')
%% load the distance bound titration and make the heatmap 
load(['output/metaboliteFPA_flexible_FPA.mat'])
targetRxns = regexprep(targetRxns,'^NewMet_','');
%% calculate correlation
dorders = 0:0.5:40;
rMat_prod = zeros(length(targetRxns),length(dorders));
FDRmat_prod = ones(length(targetRxns),length(dorders));
pMat_prod = ones(length(targetRxns),length(dorders));
CVmat_prod = ones(length(targetRxns),length(dorders));

rMat_cons = zeros(length(targetRxns),length(dorders));
FDRmat_cons = ones(length(targetRxns),length(dorders));
pMat_cons = ones(length(targetRxns),length(dorders));
CVmat_cons = ones(length(targetRxns),length(dorders));
for nn = 1: length(dorders)
    %%
    FP = FP_collection_2{1}{nn};
    relFP_prod = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
    relFP_cons = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
    for i = 1:size(FP,1)
        for j = 1:(size(FP,2)-1)
            relFP_prod(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            relFP_cons(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
        end
    end
    rmInd1 = all(isnan(relFP_prod),2) | all(relFP_prod == relFP_prod(:,1),2);
    rmInd2 = all(isnan(relFP_cons),2) | all(relFP_cons == relFP_cons(:,1),2);
    relFP_prod = relFP_prod(~rmInd1,:);
    valid_mets_prod = regexprep(targetRxns(~rmInd1),'^NewMet_','');
    relFP_cons = relFP_cons(~rmInd2,:);
    valid_mets_cons= regexprep(targetRxns(~rmInd2),'^NewMet_','');



    %Computing the correlation - production
    r=[];
    p_r=[];
    deltaminmax = [];
    testedMet = {};
    for j = 1:length(metLabel)
        metMeasure = metMat(j,:);
        if any(strcmp(valid_mets_prod,metLabel{j}))
            testedMet(end+1) = metLabel(j);
            [r(end+1),p_r(end+1)] = corr(relFP_prod(strcmp(valid_mets_prod,metLabel{j}),:)',metMeasure','type','Pearson');
            deltaminmax(end+1) = max(relFP_prod(strcmp(valid_mets_prod,metLabel{j}),:)) - min(relFP_prod(strcmp(valid_mets_prod,metLabel{j}),:));
        else
            % the FPA yeilds uniform rFP, so assumue all zeros
            testedMet(end+1) = metLabel(j);
            r(end+1) = 0;
            p_r(end+1) = 1;
            deltaminmax(end+1) = 0;
        end
    end
    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fprintf('%d mets give significant positive correlation (producing)\n',sum(r(fdr_r<0.05)>0));
    met1 = testedMet(fdr_r < 0.05 & r > 0);
    [A B] = ismember(testedMet,targetRxns);
    rMat_prod(B(A),nn) = r;
    FDRmat_prod(B(A),nn) = fdr_r;
    CVmat_prod(B(A),nn) = deltaminmax;
    pMat_prod(B(A),nn) = p_r;

    %Computing the correlation - consumption
    r=[];
    p_r=[];
    deltaminmax = [];
    testedMet = {};
    for j = 1:length(metLabel)
        metMeasure = metMat(j,:);
        if any(strcmp(valid_mets_cons,metLabel{j}))
            testedMet(end+1) = metLabel(j);
            [r(end+1),p_r(end+1)] = corr(relFP_cons(strcmp(valid_mets_cons,metLabel{j}),:)',metMeasure','type','Pearson');
            deltaminmax(end+1) = max(relFP_cons(strcmp(valid_mets_cons,metLabel{j}),:)) - min(relFP_cons(strcmp(valid_mets_cons,metLabel{j}),:));
        else
            % the FPA yeilds uniform rFP, so assumue all zeros
            testedMet(end+1) = metLabel(j);
            r(end+1) = 0;
            p_r(end+1) = 1;
            deltaminmax(end+1) = 0;
        end
    end
    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fprintf('%d mets give significant positive correlation (consuming)\n',sum(r(fdr_r<0.05)>0));
    [A B] = ismember(testedMet,targetRxns);
    rMat_cons(B(A),nn) = r;
    FDRmat_cons(B(A),nn) = fdr_r;
    CVmat_cons(B(A),nn) = deltaminmax;
    pMat_cons(B(A),nn) = p_r;
    met2 = testedMet(fdr_r < 0.05 & r > 0);

    fprintf('%d mets give significant positive correlation (union)\n',length(union(met1, met2)));
end

%% merge the results and write out tables
% default FPA
[A B] = ismember(targetRxns,localFPA_testedMets_prod);
ctrVp2 = ones(length(targetRxns),1);
ctrVr2 = zeros(length(targetRxns),1);
ctrVdeltaminmax2 = zeros(length(targetRxns),1);
ctrVp2(A) = localFPA_fdr_prod(B(A));
ctrVr2(A) = localFPA_r_prod(B(A));
ctrVdeltaminmax2(A) = localFPA_deltaminmax_prod(B(A));
% merged
pMat_producing = [ctrVp2, FDRmat_prod];
rMat_producing = [ctrVr2, rMat_prod];
deltaminmaxmat_producing = [ctrVdeltaminmax2, CVmat_prod];
t1 = array2table(rMat_producing);    
t1.Properties.RowNames = targetRxns;

% default FPA
[A B] = ismember(targetRxns,localFPA_testedMets_cons);
ctrVp2 = ones(length(targetRxns),1);
ctrVr2 = zeros(length(targetRxns),1);
ctrVdeltaminmax2 = zeros(length(targetRxns),1);
ctrVp2(A) = localFPA_fdr_cons(B(A));
ctrVr2(A) = localFPA_r_cons(B(A));
ctrVdeltaminmax2(A) = localFPA_deltaminmax_cons(B(A));
% merged
pMat_consuming = [ctrVp2, FDRmat_cons];
rMat_consuming = [ctrVr2, rMat_cons];
deltaminmaxmat_consuming = [ctrVdeltaminmax2, CVmat_cons];
t2 = array2table(rMat_consuming);    
t2.Properties.RowNames = targetRxns;
writetable(t1,'output/PCC_titration_producingPotential.csv','WriteRowNames',1);
writetable(t2,'output/PCC_titration_consumingPotential.csv','WriteRowNames',1);

%% plot the histogram of max PCC (optimal-boundary integration PCC) -- producing
r_max = [];
r_max_sig = [];
for i = 1:size(rMat_prod,1)
    if all(rMat_prod(i,:) == 0)
        r_max(i) = NaN;
        r_max_sig(i) = NaN;
    else
        r_pass = rMat_prod(i,CVmat_prod(i,:) > 0.2); % we only consider valid predictions (range > 0.2)
        fdr_pass = FDRmat_prod(i,CVmat_prod(i,:) > 0.2);
        if(length(r_pass)> 0)
            r_max(i) = max(r_pass);
            r_max_sig(i) = min(fdr_pass(r_pass == max(r_pass)));
        else
            r_max(i) = max(rMat_prod(i,:));
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
sum(r_max_sig<0.05 & r_max >0)
r_max_prod = r_max;
r_max_sig_prod = r_max_sig;
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
plt.export('figures/optimal_boundary_FPA_metabolite_correlation_producingPotential.pdf');

%% plot the histogram of max PCC (optimal-boundary integration PCC) -- consuming
r_max = [];
r_max_sig = [];
for i = 1:size(rMat_cons,1)
    if all(rMat_cons(i,:) == 0)
        r_max(i) = NaN;
        r_max_sig(i) = NaN;
    else
        r_pass = rMat_cons(i,CVmat_cons(i,:) > 0.2); % we only consider valid predictions (range > 0.2)
        fdr_pass = FDRmat_cons(i,CVmat_cons(i,:) > 0.2);
        if(length(r_pass)> 0)
            r_max(i) = max(r_pass);
            r_max_sig(i) = min(fdr_pass(r_pass == max(r_pass)));
        else
            r_max(i) = max(rMat_cons(i,:));
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
sum(r_max_sig<0.05 & r_max >0)
r_max_cons = r_max;
r_max_sig_cons = r_max_sig;
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
plt.export('figures/optimal_boundary_FPA_metabolite_correlation_consumingPotential.pdf');

%% merged
r_max = max(r_max_prod, r_max_cons);
r_max_sig = min(r_max_sig_prod, r_max_sig_cons);
figure;
hold on
histogram(r_max,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(r_max(r_max_sig<0.05 & r_max >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
hold off
sum(r_max_sig<0.05 & r_max >0)
r_max_cons = r_max;
r_max_sig_cons = r_max_sig;
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
plt.export('figures/optimal_boundary_FPA_metabolite_correlation_bestcorr.pdf');
%% filter the PCC matrix and only keep the predicted reactions 
% keep = any(pMat_valid < 0.05 & rMat_valid > 0 & deltaminmaxmat_valid > 0.2,2);
% 
% rMat_valid = rMat_valid(keep,:);
% pMat_valid = pMat_valid(keep,:);
% targetRxns_valid = targetRxns(keep);
% deltaminmaxmat_valid = deltaminmaxmat_valid(keep,:);
% 
% % row-wise normalization
% rMat_valid_normalized = rMat_valid ./ max(rMat_valid,[],2); % relative r 
% rMat_valid_normalized_ori = rMat_valid_normalized;
% 
% % label the significantly correlated predictions 
% sigMat = pMat_valid < 0.05 & rMat_valid > 0 & deltaminmaxmat_valid > 0.2;

% p-value is used here, we should 
t = struct();
t.rxn1 = targetRxns(any(FDRmat_cons < 0.05 & rMat_cons > 0 & CVmat_cons > 0.2,2));
t = struct2table(t);
writetable(t,'output/consuming_sig_rxns.csv');
t = struct();
t.rxn1 = targetRxns(any(FDRmat_prod < 0.05 & rMat_prod > 0 & CVmat_prod > 0.2,2));
t = struct2table(t);
writetable(t,'output/producing_sig_rxns.csv');

%%
% %% write out the PCC matrix whose distance label was converted to real distance 
% % bin the PCC according to converted distance 
% dorders = 0:40;
% rMat_valid_normalized = rMat_valid_normalized_ori;
% load('output/cvtDistMat.mat');
% cvtDistMat = cvtDistMat(keep,:);
% rMat_valid_normalized_binned = nan(size(rMat_valid_normalized,1),length(dorders)+2);
% sigMat_binned = nan(size(rMat_valid_normalized,1),length(dorders)+2);
% for i = 1:(length(dorders)-1)
%     lb = dorders(i);
%     ub = dorders(i+1);
%     for j = 1:size(rMat_valid_normalized_binned,1)
%         pass = find(cvtDistMat(j,:) >= lb & cvtDistMat(j,:) < ub)+2;
%         if any(pass)
%             rMat_valid_normalized_binned(j,i+2) = max(rMat_valid_normalized(j,pass));
%             sigMat_binned(j,i+2) = any(pMat_valid(j,pass) < 0.05 & rMat_valid(j,pass) > 0 & deltaminmaxmat_valid(j,pass) > 0.2);
%         end
%     end
% end
% rMat_valid_normalized_binned(:,1:2) = rMat_valid_normalized(:,1:2);
% sigMat_binned(:,1:2) = sigMat(:,1:2);
% % save data 
% t = array2table(rMat_valid_normalized_binned);
% t.Properties.RowNames = targetRxns_valid;
% writetable(t,'output/relCorr_heatmapTbl_realDist.csv','WriteRowNames',1);
% 
% IDs = [{'base 2 - boundary 6','expression only'},strsplit(num2str(dorders))];
% boundaries =  cell2table(IDs');
% writetable(boundaries,'output/heatmapTbl_boundaries_realDist.csv');
% 
% t = array2table(sigMat_binned);
% t.Properties.RowNames = targetRxns_valid;
% writetable(t,'output/heatmapTbl_sigLabel_realDist.csv','WriteRowNames',1);