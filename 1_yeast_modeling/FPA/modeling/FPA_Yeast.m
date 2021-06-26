%% notes
% flux is in  (moles / hr / mL cells); could be further normalized to
% mmole/hr/gDW by chemostat info: gDCW/ml. Should try both and see
% correlation change
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/
initCobraToolbox
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

% decouple energy from growth
% model.S(ismember(model.mets,{'s_0434[c]','s_0803[c]','s_0394[c]','s_0794[c]','s_1322[c]'}), strcmp('r_4041',model.rxns)) = 0; 

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
distMat = distMat_min;
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
fluxMat_normalized = fluxMat;

% dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
% factor = repmat(dwTbl.gDCW_mL',size(fluxMat,1),1);
% fluxMat_normalized = fluxMat_normalized * 1000 ./ factor; %mmoles/hr/gDW

GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
%  dilutionFactor = repmat([0.05 0.1 0.16 0.22 0.30],size(fluxMat,1),5);
%  fluxMat_double_normalized = fluxMat_normalized ./ dilutionFactor;
labels2 = regexprep(conditions,'_','-');
%% load merged levels
load('output/normalizedLevels_partialExcluded.mat');
%% precalc penalty 
% penalty_merged = ones(length(model.rxns),size(fluxMat,2)+1);
% [A B] = ismember(model.rxns,valid_rxns_merged);
% penalty_merged(A,1:(end-1)) = ones(size(normalizedLevel_merged(B(A),:),1),size(normalizedLevel_merged(B(A),:),2)) ./ normalizedLevel_merged(B(A),:);

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
%% original FPA
% setup some basic parameters for FPA
n = [0 1 1.5000 2 3 4 5 6 8 10]; % note: each will take 10GB mem! dont titrate all "n" and target rxns together!
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
targetRxns = rxnLabel;
% targetRxns = {'r_4041'};
%The FPA is designed with parfor loops for better speed, so we first initial the parpool
parpool(4)
[FP_collection,FP_solutions_collection] = FPA(model,targetRxns,{},distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blocklist, {},true,penalty_pro);
% [FP_collection,FP_solutions_collection] = FPA(model,targetRxns,{},distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blocklist, {},true,penalty_pro_raw);
%% FPA2
% exponential filtering could strongly filter out the influence from distal
% distance expression; and it could achive more binary filtering (lower
% than cutoff ==> strongly used; higher than distance cutoff ==> strongly
% filtered)

% by tuning the base of the exponential formula, we can define how binary
% the transition should be 
parpool(2)
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
%% check one reaction a time
targetRxn = {'r_0485'};
n2 = 18.5;
base = 1000;
[FP_collection,FP_solutions_collection] = FPA2(model,targetRxn,{},distMat,labels,n2, manualPenalty,{},max(distMat(~isinf(distMat))),blocklist, {},true,penalty_pro,base);
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

% rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
% relFP = relFP(~rmInd,:);
valid_rxns = targetRxn;%(~rmInd);
%% visualize correlation 
myLevel = relFP;
myFluxLevel = abs(fluxMat_normalized(strcmp(rxnLabel,targetRxn),:));

figure(1);
lm = fitlm(myLevel,myFluxLevel);
plot(lm);
text(myLevel,myFluxLevel,labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')
[corrR,corrP] = corr(myFluxLevel',myLevel')
%% check flux distribnution
sample = 'C0_05';
myrxn = targetRxn;
%dOrderInd = 2;
%solution = FP_solutions_collection{dOrderInd}{strcmp(myrxn,valid_rxns),strcmp(sample,conditions)}{1};
if mean(fluxMat(i,:)) > 0 
    solution = FP_solutions_collection{strcmp(myrxn,targetRxn),strcmp(sample,conditions)}{1};
else
    solution = FP_solutions_collection{strcmp(myrxn,targetRxn),strcmp(sample,conditions)}{2};
end

mytbl = listFPtbl(model,solution);
%% track flux distribution 
model_irre = convertToIrreversible(model);
%%
mytbl = listRxn(model_irre, solution.full,'s_1556[c]');
%% consistency of expression in a pathway
rxnID = 'r_0892';
printGPRForRxns(model,rxnID);
if any(strcmp(rxnID,valid_rxns_pro_perPro)) && any(strcmp(targetRxn,valid_rxns_pro_perPro))
    myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID),:);
    myLevel2 = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,targetRxn),:);
    figure(5)
    plot(myLevel,myLevel2,'.')
    text(myLevel,myLevel2,labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')
    xlabel('second rxn')
    ylabel('target rxn')
    [corrR,corrP] = corr(myLevel2',myLevel')
end
if any(strcmp(rxnID,valid_rxns_pro_perPro)) 
    myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID),:);
    myFluxLevel = abs(fluxMat_normalized(strcmp(rxnLabel,targetRxn),:));
    figure(7)
    plot(myLevel,myFluxLevel,'.')
    text(myLevel,myFluxLevel,labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')
    xlabel('second rxn')
    ylabel('flux of target rxn')
    [corrR,corrP] = corr(myFluxLevel',myLevel')
else
    fprintf('not measured!\n')
end
%% check the expression for a rxn of interest
rxnID = targetRxn;
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID),:);
myFluxLevel = abs(fluxMat_normalized(strcmp(rxnLabel,rxnID),:));

figure(2);
lm = fitlm(myLevel,myFluxLevel);
plot(lm);
text(myLevel,myFluxLevel,labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')
[corrR,corrP] = corr(myFluxLevel',myLevel')

figure(3);
X = categorical(labels2);
X = reordercats(X,labels2);
Y = myLevel;
bar(X,Y)
title(['expression of ',regexprep(rxnID,'_','')])

figure(4)
X = categorical(labels2);
X = reordercats(X,labels2);
Y = myFluxLevel;
bar(X,Y)
title(['flux of ',regexprep(rxnID,'_','')])
%% check GR
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
myrxn = {'r_2111'};
fluxMeasure = GRrate.DR_Actual';
expression = relFP(1,:);
fit = fitlm(expression,abs(fluxMeasure));
figure(1);
plot(fit)
xlabel('rFP');
ylabel('flux measurement');
[a,b] = corr(expression',abs(fluxMeasure)','type','Spearman')
labels2 = regexprep(conditions,'_','-');
text(expression,abs(fluxMeasure),labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')
%% correlation
dorders = n;
rhoMat = zeros(length(targetRxns),length(dorders));
N_sigCorr = zeros(length(dorders),1);
N_sigCorr_highCV = zeros(length(dorders),1);
perc_sigCorr_in_highCV = zeros(length(dorders),1);
for nn = 1: length(dorders)
%%
% nn = 10;
dorders(nn)
FP = FP_collection{nn};
relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
for i = 1:size(FP,1)
    for j = 1:(size(FP,2)-1)
        if  mean(fluxMat(i,:)) > 0 
            if ~isnan(FP{i,j}(1))
                relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            elseif ~isnan(FP{i,j}(2))
                relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
            else
                error('check');
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
rho=[];
p=[];
cv = [];
testedRxn = {};
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        [rho(end+1),p(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        cv(end+1) = std(relFP(strcmp(valid_rxns,rxnLabel{j}),:))/mean(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end
sum(rho(p<0.05)>0)
sum(rho(p<0.05)>0)/length(testedRxn)
histogram(rho,'BinLimits',[-1,1],'NumBins',20)
[A B] = ismember(testedRxn,targetRxns);
rhoMat(B(A),nn) = rho;
N_sigCorr(nn) = sum(rho(p<0.05)>0);
N_sigCorr_highCV(nn) = sum(rho(p<0.05 & cv > 0.1)>0);
perc_sigCorr_in_highCV(nn) = sum(rho(p<0.05 & cv > 0.1)>0) / sum(cv > 0.1);
end
%% analyze some overall metric
figure(1)
plot(dorders,N_sigCorr,'.-')
xlabel('distance order');
ylabel('number of sig corr rxns');
figure(2)
plot(dorders,N_sigCorr_highCV,'.-')
xlabel('distance order');
ylabel('number of sig-corr & rFP-CV > 0.1 rxns ');
figure(3)
plot(dorders,perc_sigCorr_in_highCV,'.-')
xlabel('distance order');
ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');
%% plot three setups together
dorders_normal = dorders/max(dorders);
N_sigCorr_normal = N_sigCorr;
N_sigCorr_highCV_normal = N_sigCorr_highCV;
perc_sigCorr_in_highCV_normal = perc_sigCorr_in_highCV;

dorders_wtdDist = dorders/max(dorders);
N_sigCorr_wtdDist = N_sigCorr;
N_sigCorr_highCV_wtdDist = N_sigCorr_highCV;
perc_sigCorr_in_highCV_wtdDist = perc_sigCorr_in_highCV;

dorders_wtdDist_expDecay = 1-dorders/max(dorders);
N_sigCorr_wtdDist_expDecay = N_sigCorr;
N_sigCorr_highCV_wtdDist_expDecay = N_sigCorr_highCV;
perc_sigCorr_in_highCV_wtdDist_expDecay = perc_sigCorr_in_highCV;
%% plot 
figure(4)
hold on
plot(dorders_normal,N_sigCorr_normal,'.-')
plot(dorders_wtdDist,N_sigCorr_wtdDist,'.-')
plot(dorders_wtdDist_expDecay,N_sigCorr_wtdDist_expDecay,'.-')
xlabel('most global --> most local (scale not comparable!)');
ylabel('number of sig corr rxns');
legend({'normal FPA','weighted distance FPA','weighted distance + exponential decay'},'Location','northwest');

figure(5)
hold on
plot(dorders_normal,N_sigCorr_highCV_normal,'.-')
plot(dorders_wtdDist,N_sigCorr_highCV_wtdDist,'.-')
plot(dorders_wtdDist_expDecay,N_sigCorr_highCV_wtdDist_expDecay,'.-')
xlabel('most global --> most local (scale not comparable!)');
ylabel('number of sig-corr & rFP-CV > 0.1 rxns ');
legend({'normal FPA','weighted distance FPA','weighted distance + exponential decay'},'Location','northwest');

figure(6)
hold on
plot(dorders_normal,perc_sigCorr_in_highCV_normal,'.-')
plot(dorders_wtdDist,perc_sigCorr_in_highCV_wtdDist,'.-')
plot(dorders_wtdDist_expDecay,perc_sigCorr_in_highCV_wtdDist_expDecay,'.-')
xlabel('most global --> most local (scale not comparable!)');
ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');
legend({'normal FPA','weighted distance FPA','weighted distance + exponential decay'},'Location','northwest');

%%
figure(2)
hold on
for i = 1:length(testedRxn)
    plot(dorders,rhoMat(i,:),'.-');
end

clustergram(rhoMat,'RowLabels',model.rxnNames(ismember(model.rxns,targetRxns)),'ColumnLabels',dorders,'Cluster','column')
%% check examples
myrxn = {'r_0473'};
specifiedExpression = {'r_0473'}
specifiedExpression = {};

%%
% r_0909 with GPR
% r_0007 without GRP
%puring pathway
% r_0908

% the new prediction: r_0473 r_1887 r_0957

fluxMeasure = fluxMat_normalized(strcmp(rxnLabel,myrxn),:);
if ~isempty(specifiedExpression)
    expression = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,specifiedExpression),:);
else
    expression = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,myrxn),:);
end
fit = fitlm(expression,abs(fluxMeasure));
figure();
plot(fit)
xlabel('expression level');
ylabel('flux measurement');
[a,b] = corr(expression',abs(fluxMeasure)','type','Spearman')
labels2 = regexprep(conditions,'_','-');
text(expression,abs(fluxMeasure),labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')

expression = relFP(strcmp(valid_rxns,myrxn),:);
fit = fitlm(expression,abs(fluxMeasure));
figure();
plot(fit)
xlabel('rFP');
ylabel('flux measurement');
[a,b] = corr(expression',abs(fluxMeasure)','type','Spearman')
text(expression,abs(fluxMeasure),labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')
%% check flux distribnution
myrxn = 'r_4041';
% mydorder = 3;
% [~,FP_solutions_collection] = FPA2(model,{myrxn},{},distMat,labels,mydorder, manualPenalty,{},max(distMat(~isinf(distMat))),blocklist, {},true,penalty_pro);
sample = 'L0_11';
solution = FP_solutions_collection_2{1,strcmp(sample,conditions)}{1};
mytbl = listFPtbl(model,solution);
%% check expression consistency 
myrxn1 = {'r_0915'};
myrxn2 = {'r_0153'};

% r_0909 with GPR
% r_0007 without GRP
%puring pathway
% r_0908

expression1 = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,myrxn1),:);
expression2 = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,myrxn2),:);

fit = fitlm(expression1,expression2);
figure();
plot(fit)
xlabel(['expression level1',myrxn1{:}]);
ylabel(['expression level2',myrxn2{:}]);
[a,b] = corr(expression',abs(fluxMeasure)','type','Spearman')
labels2 = regexprep(conditions,'_','-');
text(expression,abs(fluxMeasure),labels2,'VerticalAlignment','bottom','HorizontalAlignment','right')

%% track flux distribution 
FP_solution =  FP_solutions_collection{nn};
model_irre = convertToIrreversible(model);
%%
mytbl = listRxn(model_irre, FP_solution{strcmp(valid_rxns,myrxn),1}{1},1);
%% annotate the newly predicted rxns
newRxns = setdiff(rxns_fpa2_corr,rxns_raw_corr);
[A B] = ismember(newRxns,testedRxn);
newRxns_rho = rho(B(A));
[ROI_annotated, ROI_enrichment] = annotateRxnSet(newRxns,model);
[A B] = ismember(ROI_annotated(:,1),newRxns);
ROI_annotated(:,5) = mat2cell(newRxns_rho(B(A)),1,ones(length(newRxns),1));

size(unique(model.rxnGeneMat(ismember(model.rxns,newRxns),:),'rows'),1)





%% the permutation-based p-values (strong control of familywise error rate)
[A B] = ismember(testedRxn,rxnLabel);
fluxMeasure = fluxMat_normalized(B(A),:)';
[A B] = ismember(testedRxn,valid_rxns);
FPAmeasure = relFP(B(A),:)';

addpath mult_comp_perm_corr/
pval = [];
for i = 1:size(fluxMeasure,2)
    [pval(i), corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(fluxMeasure(:,i),FPAmeasure(:,i),100000,1,0.05,'rank');
end

[pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(fluxMeasure,FPAmeasure,100000,1,0.05,'rank');








%% find jumps and analyze the correlation with enrichment of significant correlations
normalizedLevel = relFP;
minFPcutoff = 0.8;
% then evaluate the merged dataset
AUCk_seq = 0:0.01:1;
accuracy = [];
passfilter_acc = [];
precision = [];
passfilter_prec = [];
for k = 1:length(AUCk_seq)
    isChange = min(normalizedLevel,[],2)<minFPcutoff;
    isjump = false(size(normalizedLevel,1),1);
    AUCk = AUCk_seq(k);
    for i = 1:size(normalizedLevel)
        sorted = sort(normalizedLevel(i,:));
        % scale to 0-1
        sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
        steps = sorted(2:end) - sorted(1:end-1);
        steps_sorted = sort(steps); 
        % reordered curve is
        reordered = [];
        reordered(1) = 0;
        for z = 2:length(sorted)
            reordered(z) = sum(steps_sorted(1:z-1));
        end
        AUC = trapz(1:25,reordered);
        if AUC < 12 * (1-AUCk)
            isjump(i) = true;
        else
            isjump(i) = false;
        end

    end
    noJumpRxns = valid_rxns(~isjump&isChange);
    
    % the flux data
    isjump = false(size(fluxMat,1),1);
    for i = 1:size(fluxMat)
        sorted = sort(fluxMat(i,:));
        % scale to 0-1
        sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
        steps = sorted(2:end) - sorted(1:end-1);
        steps_sorted = sort(steps); 
        % reordered curve is
        reordered = [];
        reordered(1) = 0;
        for z = 2:length(sorted)
            reordered(z) = sum(steps_sorted(1:z-1));
        end
        AUC = trapz(1:25,reordered);
        if AUC < 12 * (1-AUCk)
            isjump(i) = true;
        else
            isjump(i) = false;
        end

    end
    % perctage significant correlation with jumpping content 
    noJumpRxns2 = rxnLabel(~isjump);
    
    % use no jumping flux to access real accuracy 
    accuracy(k) = sum(p(ismember(testedRxn,noJumpRxns2))<0.05 & rho(ismember(testedRxn,noJumpRxns2)) > 0)/sum(ismember(testedRxn,noJumpRxns2));
    passfilter_acc(k) = sum(ismember(testedRxn,noJumpRxns2)) / length(testedRxn);

    % use double no jumping to access recall
    precision(k) = sum(p(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))<0.05 & rho(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))>0)/...
    sum(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)));
    passfilter_prec(k) = sum(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2))) / length(testedRxn);
end
figure;
hold on
plot(AUCk_seq,accuracy);
plot(AUCk_seq,passfilter_acc);
plot(AUCk_seq,precision);
plot(AUCk_seq,passfilter_prec);
xlim([0.3,1])
legend({'accuracy','passfilter%(acc)','precision','passfilter%(prec)'});
% we pick 0.61 as best 
% this analysis reveals that the dosage-responsively regualated enzymes are
% good proxy for flux, or the jumppings are insignificant because of lack
% of power (see below)
%% visualize seperation and correlation of expression
isChange = min(normalizedLevel,[],2)<0.8;
isjump = false(size(normalizedLevel,1),1);
AUCk = 0.59;
for i = 1:size(normalizedLevel)
    sorted = sort(normalizedLevel(i,:));
    % scale to 0-1
    sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
    steps = sorted(2:end) - sorted(1:end-1);
    steps_sorted = sort(steps); 
    % reordered curve is
    reordered = [];
    reordered(1) = 0;
    for z = 2:length(sorted)
        reordered(z) = sum(steps_sorted(1:z-1));
    end
    AUC = trapz(1:25,reordered);
    if AUC < 12 * (1-AUCk)
        isjump(i) = true;
    else
        isjump(i) = false;
    end

end
% perctage significant correlation with jumpping content 
jumpRxns = valid_rxns(isjump&isChange);
noJumpRxns = valid_rxns(~isjump&isChange);
%nochangerxn = valid_rxns(~isChange);
%sum(p(ismember(testedRxn,nochangerxn))<0.05)/sum(ismember(testedRxn,nochangerxn))
%sum(p(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))<0.05)/sum(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))
%sum(p(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))<0.05) / sum(p<0.05)
jumpping = normalizedLevel(isjump&isChange,:);
figure;
hold on 
% scale to 0-1
jumpping = (jumpping - min(jumpping,[],2)) ./ repmat((max(jumpping,[],2) - min(jumpping,[],2)),1,25);
for i = 1:size(jumpping,1)
    plot(sort(jumpping(i,:)),'r');
end
nojumpping = normalizedLevel(~isjump&isChange,:);
nojumpping = (nojumpping - min(nojumpping,[],2)) ./ repmat((max(nojumpping,[],2) - min(nojumpping,[],2)),1,25);
figure;
hold on 
for i = 1:size(nojumpping,1)
    plot(sort(nojumpping(i,:)),'g');
end
% histogram(rho(ismember(testedRxn,noJumpRxns)))


%%  seperate out the jummping fluxes (for precision)
isjump = false(size(fluxMat,1),1);
for i = 1:size(fluxMat)
    sorted = sort(fluxMat(i,:));
    % scale to 0-1
    sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
    steps = sorted(2:end) - sorted(1:end-1);
    steps_sorted = sort(steps); 
    % reordered curve is
    reordered = [];
    reordered(1) = 0;
    for z = 2:length(sorted)
        reordered(z) = sum(steps_sorted(1:z-1));
    end
    AUC = trapz(1:25,reordered);
    if AUC < 12 * (1-AUCk)
        isjump(i) = true;
    else
        isjump(i) = false;
    end

end
% perctage significant correlation with jumpping content 
jumpRxns2 = rxnLabel(isjump);
noJumpRxns2 = rxnLabel(~isjump);

jumpping = fluxMat(isjump,:);
figure;
hold on 
% scale to 0-1
jumpping = (jumpping - min(jumpping,[],2)) ./ repmat((max(jumpping,[],2) - min(jumpping,[],2)),1,25);
for i = 1:size(jumpping,1)
    plot(sort(jumpping(i,:)));
end
nojumpping = fluxMat(~isjump,:);
nojumpping = (nojumpping - min(nojumpping,[],2)) ./ repmat((max(nojumpping,[],2) - min(nojumpping,[],2)),1,25);
figure;
hold on 
for i = 1:size(nojumpping,1)
    plot(sort(nojumpping(i,:)));
end

% use no jumping flux to access real accuracy == now we only use positive
% correlation as accessing prediction accuracy
sum(p(ismember(testedRxn,noJumpRxns2))<0.05 & rho(ismember(testedRxn,noJumpRxns2)) > 0)/sum(ismember(testedRxn,noJumpRxns2))
%sum(p(ismember(testedRxn,noJumpRxns2))<0.05) / sum(p<0.05)
sum(ismember(testedRxn,noJumpRxns2)) / length(p)
figure;
histogram(rho(ismember(testedRxn,noJumpRxns2)))

% use double no jumping to access recall
sum(p(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))<0.05 & rho(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))>0)/...
    sum(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))
sum(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2))) / length(p)
%sum(p(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2)))<0.05) / sum(p<0.05)
figure;
histogram(rho(ismember(testedRxn,intersect(noJumpRxns,noJumpRxns2))))



% find a linear pathway and check what is going wrong with FPA
% also individually analyze the 56 SIMMER reactiosn


















%% precision-recall curve
master_expression =master_expression_pro;
nn = 4;
n(nn)
rho_f=nan(length(rxnLabel),1);
p_f=nan(length(rxnLabel),1);
rho_r=nan(length(rxnLabel),1);
p_r=nan(length(rxnLabel),1);
FP = FP_collection{nn};
relFP_f = nan(size(FP,1),length(master_expression));%flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));%flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
cutoffs= 0:0.02:1;
binCallK = 0.3; %up or down 25%
z_fluxMat = normalize(fluxMat,2,'zscore');
precision = nan(length(rxnLabel),length(cutoffs));
recall = nan(length(rxnLabel),length(cutoffs));
binaryObs = {};
for z = 1:length(rxnLabel)
    if  mean(fluxMat(z,:)) > 0 && ~isnan(relFP_f(z,1))
        [precision(z,:),recall(z,:),binaryObs{z}] = prCurve(relFP_f(z,:),z_fluxMat(z,:),cutoffs,binCallK);
    else
        [precision(z,:),recall(z,:),binaryObs{z}] = prCurve(relFP_r(z,:),z_fluxMat(z,:),cutoffs,binCallK);
    end
end
figure;
hold on
for ii = 1:length(rxnLabel)
    plot(recall(ii,:),precision(ii,:))
    plot(sum(binaryObs{ii}==1) / sum(binaryObs{ii}==1 | binaryObs{ii}==-1),sum(binaryObs{ii}==1) / length(binaryObs{ii}),'r.','MarkerSize',10)
    plot(sum(binaryObs{ii}==-1) / sum(binaryObs{ii}==1 | binaryObs{ii}==-1),sum(binaryObs{ii}==-1) / length(binaryObs{ii}),'r.','MarkerSize',10)
end
%refline(1,0)
hold off
% for z = 1:length(rxnLabel)
%     baseline(z) = baselinePrecision(abs(fluxMat(z,:)),binCallK*abs(fluxMat(z,:)));
% end
% hold on
% for ii = 1:length(rxnLabel)
%  refline(0,baseline(ii))
% end
% internal only
figure;
hold on
indexes = find(internalRxnInd);
for ii = 1:length(indexes)
 plot(recall(indexes(ii),:),precision(indexes(ii),:))
end
refline(1,0)
% external only
figure;
hold on
indexes = find(~internalRxnInd);
for ii = 1:length(indexes)
 plot(recall(indexes(ii),:),precision(indexes(ii),:))
end
refline(1,0)
%% overall precision recall
for z = 1:length(rxnLabel)
    if  mean(fluxMat(z,:)) > 0 && ~isnan(relFP_f(z,1))
        relFP(z,:) = relFP_f(z,:);
    else
        relFP(z,:) = relFP_r(z,:);
    end
end
% if normalize by diluition rate
relFP_normalized = relFP;
group1 = [1,6,11,16,21];
group2 = [2,7,12,17,22];
group3 = [3,8,13,18,23];
group4 = [4,9,14,19,24];
group5 = [5,10,15,20,25];
relFP_normalized(:,group1) = relFP_normalized(:,group1) ./ repmat(max(abs(relFP_normalized(:,group1)),[],2),1,5);
relFP_normalized(:,group2) = relFP_normalized(:,group2) ./ repmat(max(abs(relFP_normalized(:,group2)),[],2),1,5);
relFP_normalized(:,group3) = relFP_normalized(:,group3) ./ repmat(max(abs(relFP_normalized(:,group3)),[],2),1,5);
relFP_normalized(:,group4) = relFP_normalized(:,group4) ./ repmat(max(abs(relFP_normalized(:,group4)),[],2),1,5);
relFP_normalized(:,group5) = relFP_normalized(:,group5) ./ repmat(max(abs(relFP_normalized(:,group5)),[],2),1,5);

%%
binCallK = 0.3;
z_fluxMat_double_normalized = normalize(fluxMat_double_normalized,2,'zscore');
[precision,recall,binaryObs] = prCurve_lump(relFP,z_fluxMat_double_normalized,cutoffs,binCallK);
figure;
plot(recall,precision)
hold on
plot(sum(binaryObs==1) / sum(binaryObs==1 | binaryObs==-1),sum(binaryObs==1) / length(binaryObs),'r.','MarkerSize',10)
plot(sum(binaryObs==-1) / sum(binaryObs==1 | binaryObs==-1),sum(binaryObs==-1) / length(binaryObs),'r.','MarkerSize',10)
%%
%% normalize by dilution rate
fluxMat_double_normalized = fluxMat_normalized;
group1 = [1,6,11,16,21];
group2 = [2,7,12,17,22];
group3 = [3,8,13,18,23];
group4 = [4,9,14,19,24];
group5 = [5,10,15,20,25];
fluxMat_double_normalized(:,group1) = fluxMat_double_normalized(:,group1) ./ repmat(max(abs(fluxMat_double_normalized(:,group1)),[],2),1,5);
fluxMat_double_normalized(:,group2) = fluxMat_double_normalized(:,group2) ./ repmat(max(abs(fluxMat_double_normalized(:,group2)),[],2),1,5);
fluxMat_double_normalized(:,group3) = fluxMat_double_normalized(:,group3) ./ repmat(max(abs(fluxMat_double_normalized(:,group3)),[],2),1,5);
fluxMat_double_normalized(:,group4) = fluxMat_double_normalized(:,group4) ./ repmat(max(abs(fluxMat_double_normalized(:,group4)),[],2),1,5);
fluxMat_double_normalized(:,group5) = fluxMat_double_normalized(:,group5) ./ repmat(max(abs(fluxMat_double_normalized(:,group5)),[],2),1,5);

%% flux distribution
model_irrev = convertToIrreversible(model); % convert to irreversible model
% to inspect the flux distribution, we provide a simple flux tracker 
mytbl = listRxn(model_irrev,FP_solutions{3,1}{1}.full,'glu_L[e]')
%% new FPA
% we only compare the first 10 cell lines for the sake of time
% master_expression = master_expression(1:10);
% setup some basic parameters for FPA
model = changeRxnBounds(model,'r_4046',1,'b'); % maintance 
% define the flux nutrient
P_total = abs(sum(fluxMat_normalized(ismember(rxnLabel,{'r_1244'}),:),1));
C_total = abs(sum(fluxMat_normalized(ismember(rxnLabel,{'r_1166'}),:),1));
N_total = abs(sum(fluxMat_normalized(ismember(rxnLabel,{'r_1115'}),:),1));
U_total = abs(sum(fluxMat_normalized(ismember(rxnLabel,{'r_1272'}),:),1));
L_total = abs(sum(fluxMat_normalized(ismember(rxnLabel,{'r_1211'}),:),1));
% make master constriants ==> give different constriant for different
% conditions is wired!
% master_const = {};% we call this variable "master_expression"
% for i = 1:length(conditions)
%     constriants = struct();
%     constriants.rxns = {'r_2005','r_1714','r_1654','r_2090','r_1899'};
%     % phosphate, glucose, ammonium, uracil, leucine
%     str = conditions{i};
%     switch str(1)
%         case 'P'
%             constriants.lb = [];
%             constriants.ub = [1000,1000,1000,1000,1000];
%     end
%     
%     master_const{i} = constriants;
% end
% we use the glucose limiting condition with minimum dilution rate to
% construct the reference 
% phosphate exchange
model.lb(strcmp(model.rxns,'r_2005')) = -P_total(strcmp(conditions,'C0_05'));
% glucose exchange
model.lb(strcmp(model.rxns,'r_1714')) = -C_total(strcmp(conditions,'C0_05'));
% ammonium exchange 
model.lb(strcmp(model.rxns,'r_1654')) = -N_total(strcmp(conditions,'C0_05'));
% uracil 
model.lb(strcmp(model.rxns,'r_2090')) = 0;
% leucine
model.lb(strcmp(model.rxns,'r_1899')) = 0;

n = [0,0.5,1.5,3,10];
changeCobraSolverParams('LP','optTol', 10e-8);
changeCobraSolverParams('LP','feasTol', 10e-8);
targetRxns = rxnLabel;
%The FPA is designed with parfor loops for better speed, so we first initial the parpool
[FP_collection,FP_solutions_collection] = FPA_minAllowance(model,targetRxns,master_expression_pro,distMat,labels,n, manualPenalty);

%% metabolite centric
 % it is like a demand reaction except for adding a new reaction
myMet = 'lac_L[c]';
% find the associated end reactions (producing the metabolite)
myRxns = model.rxns(any(model.S(cellfun(@(x) ~isempty(regexp(x,['^',myMet(1:end-3),'_.\[.\]$'],'once')),model.mets),:),1)); %all rxns use these mets
myRxns = myRxns(cellfun(@(x) ~isempty(regexp(x,'^(sink_|EX_)','once')),myRxns));
% block these reactions and optimize 
% only block the reactions for corresponding tissue
model_tmp = model;
% add the new demand 
model_tmp = addReaction(model_tmp,['DMN',myMet],'reactionFormula',[myMet,' ->'],'geneRule', '','printLevel',0);
targetrxn_fullName = ['DMN',myMet];
% the distance of this new demand rxn needs to be calculated; Note
% that the distance is actually 1+min(distance(any rxn that
% contains this met)
[distMat_raw_2,labels_2] = updateDistance(targetrxn_fullName,model_tmp,distMat_raw,labels,max(distMat(~isinf(distMat))));
% take the min and block intestine
distMat_2 = distMat_raw_2;
for ii = 1:size(distMat_2,1)
    for jj = 1:size(distMat_2,2)
        distMat_2(ii,jj) = min([distMat_raw_2(ii,jj),distMat_raw_2(jj,ii)]);
    end
end
model_tmp.ub(ismember(model_tmp.rxns,myRxns)) = 0;
model_tmp.lb(ismember(model_tmp.rxns,myRxns)) = 0;
[FP,FP_solutions] = FPA_minAllowance(model_tmp,{targetrxn_fullName},master_expression,distMat_2,labels_2,n, manualPenalty);
%% 5. make relative flux potential
relFP_f = nan(size(FP,1),length(master_expression));%flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));%flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
%% 6. compare with measured lactate production
qryRxn = targetRxns{1}
lookuptbl = readtable('input/palssonIDtbl.xlsx');
[A B] = ismember(conditions,lookuptbl.expID);
polID = lookuptbl.palssonID(B(A));
validPredInd = A;
[A B] = ismember(polID,cellLabel);
meanFlux = (flux_lb + flux_ub)./2;
Lact = meanFlux(strcmp(rxnLabel,qryRxn),B(A));
fit = fitlm(Lact,relFP_f(strcmp(targetRxns,qryRxn),validPredInd));
plot(fit)


%%
x = Lact;
y = relFP_f(strcmp(targetRxns,qryRxn),validPredInd);
removeInd = (x <=0.9 & y>=0.65) | (x >= 1.5 & y<=0.6);
fit = fitlm(x(~removeInd),y(~removeInd));
plot(fit)
fit
failToPred = polID(removeInd);
name2 = conditions(ismember(conditions,lookuptbl.expID));
name2(removeInd)'

%% show correlation with gene expression
myGene = {'ENSG00000163738','ENSG00000065911'};
qryRxn = 'HMR_9042';
geneVec = expTbl(ismember(expTbl.Ensembl,myGene),:);
lookuptbl = readtable('input/palssonIDtbl.xlsx');
[A B] = ismember(conditions,lookuptbl.expID);
validPredInd = A;
polID = lookuptbl.palssonID(B(A));
[A B] = ismember(polID,cellLabel);
meanFlux = (flux_lb + flux_ub)./2;
Flux = abs(meanFlux(strcmp(rxnLabel,qryRxn),B(A)));
[A B] = ismember(conditions(validPredInd),geneVec.Properties.VariableNames);
fit = fitlm(sum(geneVec{:,B(A)},1),Flux);
plot(fit)
%% correlation for all
Data = struct();
Data.gene_name = expTbl.Ensembl;
Data.GE = expTbl{:,4:end};
[A B] = ismember(conditions,lookuptbl.expID);
validPredInd = A;
Data.GE(:,~validPredInd)=[];
Data.GR = Lact;
%% predtable for some cell type
sel = [1:5,12:18,25:34,51:52];
fit = fitlm(exp(sel),relFP_f(sel));
fit
plot(fit)

%% 6. compare the growth rate
GR = readtable('GrowthRate.xlsx');
exp = [];
predSampleName = expTbl.Properties.VariableNames(3:12);
for i = 1:length(predSampleName)
    exp(i) = GR.GR(strcmp(GR.CellLine_ID_ORI,predSampleName{i}));
end
fit = fitlm(exp,relFP_f(2,:));
plot(fit)
%% old codes
%% 1. inspect the flux distribution for reported FP values
% the flux distribution is reported for the irreversible model, in the
% FP_solutions. 
% to get irreversible model
model_irrev = convertToIrreversible(model); % convert to irreversible model
% to inspect the flux distribution, we provide a simple flux tracker 
mytbl = listRxn(model_irrev,FP_solutions{3,1}{1}.full,'glu_L[e]')
%%
model_unconst_irrev = convertToIrreversible(model_unconst); % convert to irreversible model
% to inspect the flux distribution, we provide a simple flux tracker 
mytbl = listRxn(model_unconst_irrev,FP_solutions{3,4}{1}.full,'glu_L[c]')
%% 
model_sNutr_irrev = convertToIrreversible(model_sNutr); % convert to irreversible model
% to inspect the flux distribution, we provide a simple flux tracker 
mytbl = listRxn(model_sNutr_irrev,FP_solutions{2,61}{1}.full,'lac_L[c]')
%% evaluate n titrate
master_expression =master_expression_pro;
nn = 4;
n(nn)
rho_f=nan(length(rxnLabel),1);
p_f=nan(length(rxnLabel),1);
rho_r=nan(length(rxnLabel),1);
p_r=nan(length(rxnLabel),1);
FP = FP_collection{nn};
relFP_f = nan(size(FP,1),length(master_expression));%flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));%flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
for z = 1:length(rxnLabel)
    if  ~isnan(relFP_f(z,1))
        [rho_f(z),p_f(z)] = corr(relFP_f(z,:)',abs(fluxMat_double_normalized(z,:))','type','Spearman');
    else
        [rho_f(z),p_f(z)] = deal(nan,nan);
    end
    if  ~isnan(relFP_r(z,1))
        [rho_r(z),p_r(z)] = corr(relFP_r(z,:)',abs(fluxMat_double_normalized(z,:))','type','Spearman');
    else
        [rho_r(z),p_r(z)] = deal(nan,nan);
    end
end
best_rho = max([rho_f,rho_r],[],2);
best_p = min([p_f,p_r],[],2);
sum(best_p<0.05)
sum(best_p<0.05) / length(rxnLabel)
sum(best_p<0.05 & best_rho>0)
sum(best_p<0.05 & best_rho>0) / length(rxnLabel)
figure(3)
histogram(best_p)
figure(4)
histogram(best_rho)
%% only compare those with direct expression values
penalty = calculatePenalty(model,master_expression_pro);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only data-derived level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns_pro = model.rxns(~rmInd);
comparableRxns_pro = intersect(rxnLabel,valid_rxns_pro);
%% compare
sum(best_p<0.05 & ismember(rxnLabel,comparableRxns_pro))
sum(best_p<0.05 & ismember(rxnLabel,comparableRxns_pro)) / length(comparableRxns_pro)
figure(3)
histogram(log10(best_p(ismember(rxnLabel,comparableRxns_pro))))
figure(4)
histogram(best_rho(ismember(rxnLabel,comparableRxns_pro)))
% no improvement at all? 
%% plot example 
figure;
fluxQry = 'r_0959';%r_0816 r_0959 r_0152 r_0727 r_0889
fluxMeasure = fluxMat_double_normalized(strcmp(rxnLabel,fluxQry),:);
FPAMeasure = relFP_f(strcmp(rxnLabel,fluxQry),:);
corr(FPAMeasure',abs(fluxMeasure)','type','Spearman')
fit = fitlm(FPAMeasure',abs(fluxMeasure)');
plot(fit)
