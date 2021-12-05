%% About
% the case study of several reactions 
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
%% 3. load the distance matrix, and other optional inputs
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
% make the raw flux matrix (in per gDW unit)
dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
fluxMat_raw = fluxMat;
factor = repmat(dwTbl.gDCW_mL',size(fluxMat_raw,1),1);
fluxMat_raw = fluxMat_raw * 1000 ./ factor; %mmoles/hr/gDW

%% load expression levels
load('output/normalizedLevels_partialExcluded.mat');
%% precalc penalty 
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

%% load the FPA data 
load(['output/Titration_relativeExp_wtdDist_expDecay.mat'])
%% correlation for 2-d titration
zz = 2;
nn = 9;
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

sigCorrRxns_FPA = testedRxn(fdr_r<0.05 & r >0  & deltaminmax > 0.2);
corrTbl_FPA = table(testedRxn',r',fdr_r');
%% check the expression for a rxn of interest
%% case 1: the proline biosynthesis pathway
% r_1887, up r_0473, down r_0957
rxnID1 = 'r_0473';
rxnID2 = 'r_1887';
rxnID3 = 'r_0957';
rxnID4 = 'r_0468';

figure;
hold on
myFlux = fluxMat_normalized(strcmp(rxnLabel,rxnID2),:);
myFlux = normalize(myFlux,'range',[0.01 1]);
[~, order0] = sort(myFlux);

myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID1),:);
X = categorical(conditions(order0));
X = reordercats(X,conditions(order0));
Y = myLevel(order0);
Y1 = normalize(Y,'range',[0.01 1]);
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID3),:);
Y = myLevel(order0);
Y2 = normalize(Y,'range',[0.01 1]);
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID4),:);
Y = myLevel(order0);
Y3 = normalize(Y,'range',[0.01 1]);
b = bar(X,[Y3',Y1',Y2']);
b(1).FaceColor = '#D95319';
b(2).FaceColor = '#0072BD';
b(3).FaceColor = '#EDB120';

% flux
plot(X,myFlux(order0),'o-','Color','#77AC30','LineWidth',3)
% flux potential 
myFP = relFP(strcmp(valid_rxns,rxnID2),:);
Y = myFP(order0);
Y = normalize(Y,'range',[0.01 1]);
plot(X,Y,'o-','Color','#4DBEEE','LineWidth',3)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.5, 4];
plt.LineWidth = [1 1 1 1.5 1.5];
plt.FontName = 'Arial';
legend({'relative expression of r_0468','relative expression of r_0473','relative expression of r_0957','relative flux of r_1887','rFP of r_1887'},'FontSize',12)
ylabel('Relative level',  'FontSize',7)
plt.Interpreter = 'none';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.LegendLoc = 'eastoutside';
plt.FontSize = 7;
plt.export('figures/FPA_showcase_r_1887_mechanism.pdf');
%% make the figure of other reactions
%% load the FPA data 
load(['output/Titration_relativeExp_oriDist_oriDecay.mat'])
rxnID1 = 'r_0565'; % also for r_0565
nn = 4;
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

figure;
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID1),:);
myFlux = fluxMat_normalized(strcmp(rxnLabel,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','.')
set(h(1),'MarkerSize',15)
set(h(2), {'color'},{'#808080'}) 
set(h(2), {'LineStyle'},{'--'}) 
set(h(3), {'visible'},{'off'}) 
set(h(4), {'visible'},{'off'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
load('output/sig_corr_rxns_rel2rel.mat');
text(0.4,max(myFlux)*0.93,['r = ',num2str(corrTbl.Var2(strcmp(corrTbl.Var1,rxnID1)),2)],'FontSize',15)
text(0.4,max(myFlux)*0.85,['FDR = ',num2str(corrTbl.Var3(strcmp(corrTbl.Var1,rxnID1)),2)],'FontSize',15)

plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Title = rxnID1;
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';

plt.export(['figures/others/flux_expression_',rxnID1,'.pdf']);

figure;
myLevel = relFP(strcmp(valid_rxns,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','.')
set(h(1),'MarkerSize',15)
set(h(2), {'color'},{'#808080'}) 
set(h(2), {'LineStyle'},{'--'}) 
set(h(3), {'visible'},{'off'}) 
set(h(4), {'visible'},{'off'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
text(0.47,max(myFlux)*0.93,['r = ',num2str(r(strcmp(testedRxn,rxnID1)),2)],'FontSize',15)
text(0.47,max(myFlux)*0.85,['FDR = ',num2str(fdr_r(strcmp(testedRxn,rxnID1)),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Title = rxnID1;
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.export(['figures/others/original_rFP_expression_',rxnID1,'.pdf']);

% third plot
load(['output/Titration_relativeExp_wtdDist_expDecay.mat']);
zz = 2;
nn = 9;
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
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & deltaminmax > 0.1)>0));


figure;
myLevel = relFP(strcmp(valid_rxns,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
xlim1 = xlim;
ylim1 = ylim;
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','.')
set(h(1),'MarkerSize',15)
set(h(2), {'color'},{'#808080'}) 
set(h(2), {'LineStyle'},{'--'}) 
set(h(3), {'visible'},{'off'}) 
set(h(4), {'visible'},{'off'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
text(0.6,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
text(0.6,max(myFlux)*0.85,['FDR = ',num2str(fdr_r(strcmp(testedRxn,rxnID1)),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Title = rxnID1;
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.XLim = xlim1;
plt.YLim = ylim1;
plt.export(['figures/others/new_rFP_expression_',rxnID1,'.pdf']);
