addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/

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
% load proteomics
% proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
% % we use the non-log level
% proTbl{:,2:end} = 2.^proTbl{:,2:end};
% % preprocess the expression table
% % to facilate the future use of the expression of many samples, we
% % re-organize it into a structure variable.
% % the FPA matrix will be in the same order as the master_expression
% conditions = proTbl.Properties.VariableNames(2:end);
% % make a new master_expression for these four conditions.
% master_expression = {};% we call this variable "master_expression"
% geneInd = ismember(proTbl.Gene, model.genes); % get the index of genes in the model
% for i = 1:length(conditions)
%     expression = struct();
%     expression.genes = proTbl.Gene(geneInd);
%     expression.value = proTbl.(conditions{i})(geneInd);
%     master_expression{i} = expression;
% end
% master_expression_pro = master_expression;
% 
% % load microarray
% rnaTbl = readtable('./../input/YeastJoshua/MicroArray/matched_knn_imputed_log2_FC_to_reference_pmid_17959824.txt');% this is the log2(FC_reference)
% % we use the non-log level
% rnaTbl{:,3:end} = 2.^rnaTbl{:,3:end};
% % preprocess the expression table
% % to facilate the future use of the expression of many samples, we
% % re-organize it into a structure variable.
% % the FPA matrix will be in the same order as the master_expression
% % make a new master_expression for these four conditions.
% master_expression = {};% we call this variable "master_expression"
% geneInd = ismember(rnaTbl.YORF, model.genes); % get the index of genes in the model
% for i = 1:length(conditions)
%     expression = struct();
%     expression.genes = rnaTbl.YORF(geneInd);
%     expression.value = rnaTbl.(conditions{i})(geneInd);
%     master_expression{i} = expression;
% end
% master_expression_rna = master_expression;

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
% normalize flux unit
% when normalize the flux by the flux of biomass production (growth rate),
% the unit of growth rate needs to be taken care of. In chemostat setting,
% steady state was defined as stable OD (see SIMMER paper), which means
% steady cell density (number, aka, volume). Therefore, the dilution rate
% is a measure of per cell flux. So, we should normalize the internal flux
% under /ml cell metric

% flux is in  (moles / hr / mL cells); no conversion is needed. 
% in fact, correlation got worse if we normalzie the flux to / gDW first!
fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
%  dilutionFactor = repmat([0.05 0.1 0.16 0.22 0.30],size(fluxMat,1),5);
%  fluxMat_double_normalized = fluxMat_normalized ./ dilutionFactor;

% make the raw flux matrix (in per gDW unit)
% flux is in  (moles / hr / mL cells); could be further normalized to
% mmole/hr/gDW by chemostat info: gDCW/ml. such that it is comparable with
% the per gDW protein content
dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
fluxMat_raw = fluxMat;
factor = repmat(dwTbl.gDCW_mL',size(fluxMat_raw,1),1);
fluxMat_raw = fluxMat_raw * 1000 ./ factor; %mmoles/hr/gDW

%% find internal reactions
% EXrelRxns = model.rxns(any(model.S(cellfun(@(x) ~isempty(regexp(x,'\[e\]$','once')),model.mets),:),1));
% transporters = intersect(rxnLabel,EXrelRxns);
% printRxnFormula_XL(model,transporters);
% internalRxnInd = ~ismember(rxnLabel,transporters);
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
%Computing the correlation between a reaction expression and measured growth rate
r=[];
p_r=[];
cv = [];
testedRxn = {};
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        cv(end+1) = std(relFP(strcmp(valid_rxns,rxnLabel{j}),:))/mean(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end
            
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & cv > 0.1)>0));

sigCorrRxns_FPA = testedRxn(fdr_r<0.05 & r >0  & cv > 0.1);
corrTbl_FPA = table(testedRxn',r',fdr_r');

%% annotate the newly predicted rxns
load('output/sig_corr_rxns_raw2raw.mat');
sigCorrRxns_raw2raw = sigCorrRxns;
corrTbl_raw2raw = corrTbl;
load('output/sig_corr_rxns_rel2rel.mat');
sigCorrRxns_rel2rel = sigCorrRxns;
corrTbl_rel2rel = corrTbl;

newRxns = setdiff(sigCorrRxns_FPA,intersect(sigCorrRxns_FPA,sigCorrRxns_rel2rel));
[ROI_annotated, ROI_enrichment] = annotateRxnSet(newRxns,model,testedRxn);

ROI_annotated(:,5) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl_rel2rel.Var1);
ROI_annotated(A,5) = mat2cell(corrTbl_rel2rel{B(A),2},ones(sum(A),1),1); % pearson r

ROI_annotated(:,6) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl_FPA.Var1);
ROI_annotated(A,6) = mat2cell(corrTbl_FPA{B(A),2},ones(sum(A),1),1); % pearson r

ROI_annotated(:,7) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl_rel2rel.Var1);
ROI_annotated(A,7) = mat2cell(corrTbl_rel2rel{B(A),3},ones(sum(A),1),1); % p

ROI_annotated(:,8) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl_FPA.Var1);
ROI_annotated(A,8) = mat2cell(corrTbl_FPA{B(A),3},ones(sum(A),1),1); % p

ROI_annotated = [{'rxnID','formula','GPR','firstSubsystem','exp2flux pearson r','exp2flux FDR','rFP2flux pearson r','rFP2flux FDR'};ROI_annotated];
ROI_enrichment = [{'SubSystems','number of significantly correlated rxns','FDR'};ROI_enrichment];

% ROI_annotated = cell2table(ROI_annotated);
size(unique(model.rxnGeneMat(ismember(model.rxns,newRxns),:),'rows'),1)

% fig = uifigure;
% uit = uitable(fig,'Data',ROI_annotated,'ColumnName',ROI_annotated.Properties.VariableNames,...
%     'RowName',ROI_annotated.Properties.RowNames);
%% check the expression for a rxn of interest
%% case 1
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
% another gene expression
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID3),:);
Y = myLevel(order0);
Y2 = normalize(Y,'range',[0.01 1]);
% b = bar(X,[Y1',Y2']);
% b(1).FaceColor = '#D95319';
% b(2).FaceColor = '#0072BD';
% the third gene expression
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
plt.BoxDim = [13, 8];
plt.LineWidth = [1.5 1.5 1.5 2.5 2.5];
plt.FontSize = 12;
plt.FontName = 'Arial';
legend({'relative expression of r_0468','relative expression of r_0473','relative expression of r_0957','relative flux of r_1887','rFP of r_1887'},'FontSize',12)
ylabel('Relative level',  'FontSize',15)
plt.Interpreter = 'none';
plt.LegendLoc = 'eastoutside';
plt.export('figures/FPA_showcase_r_1887_mechanism.tiff');

% plot(X, normalize(mean([Y1;Y2]),'range',[0.01 1]),'go-')
%% plot the correlation of FPA-flux and exp-flux for rxn1, rxn2 and rxn3
figure;
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID4),:);
myFlux = fluxMat_normalized(strcmp(rxnLabel,rxnID4),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
text(0.3,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/FPA_showcase_r_1887_flux_expression_r_0468.tiff']);
figure;
myLevel = relFP(strcmp(valid_rxns,rxnID4),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Flux Potential');
ylabel('Relative Flux (absolute value)');
text(0.3,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/FPA_showcase_r_1887_flux_FP_r_0468.tiff']);



figure;
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID1),:);
myFlux = fluxMat_normalized(strcmp(rxnLabel,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
text(0.3,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/FPA_showcase_r_1887_flux_expression_r_0473.tiff']);
figure;
myLevel = relFP(strcmp(valid_rxns,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Flux Potential');
ylabel('Relative Flux (absolute value)');
text(0.3,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/FPA_showcase_r_1887_flux_FP_r_0473.tiff']);

figure;
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID3),:);
myFlux = fluxMat_normalized(strcmp(rxnLabel,rxnID3),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
text(0.4,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/FPA_showcase_r_1887_flux_expression_r_0957.tiff']);
figure;
myLevel = relFP(strcmp(valid_rxns,rxnID3),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Flux Potential');
ylabel('Relative Flux (absolute value)');
text(0.3,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/FPA_showcase_r_1887_flux_FP_r_0957.tiff']);

%% oratate transporter could be another case

%% show case of systems level predictions: glucose transporter


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figure of any reaction
%% load the FPA data 
load(['output/Titration_relativeExp_oriDist_oriDecay.mat'])
rxnID1 = 'r_0485';
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
%Computing the correlation between a reaction expression and measured growth rate
r=[];
p_r=[];
cv = [];
testedRxn = {};
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        cv(end+1) = std(relFP(strcmp(valid_rxns,rxnLabel{j}),:))/mean(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end
            
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & cv > 0.1)>0));

figure;
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID1),:);
myFlux = fluxMat_normalized(strcmp(rxnLabel,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
load('output/sig_corr_rxns_rel2rel.mat');
text(0.4,max(myFlux)*0.93,['r = ',num2str(corrTbl.Var2(strcmp(corrTbl.Var1,rxnID1)),2)],'FontSize',15)
text(0.4,max(myFlux)*0.85,['FDR = ',num2str(corrTbl.Var3(strcmp(corrTbl.Var1,rxnID1)),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.Title = rxnID1;
plt.Interpreter = 'None';
plt.export(['figures/others/flux_expression_',rxnID1,'.tiff']);

figure;
myLevel = relFP(strcmp(valid_rxns,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
text(0.47,max(myFlux)*0.93,['r = ',num2str(r(strcmp(testedRxn,rxnID1)),2)],'FontSize',15)
text(0.47,max(myFlux)*0.85,['FDR = ',num2str(fdr_r(strcmp(testedRxn,rxnID1)),2)],'FontSize',15)
xlabel('Relative Flux Potential');
ylabel('Relative Flux (absolute value)');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.Title = rxnID1;
plt.Interpreter = 'None';
plt.export(['figures/others/original_rFP_expression_',rxnID1,'.tiff']);

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
%Computing the correlation between a reaction expression and measured growth rate
r=[];
p_r=[];
cv = [];
testedRxn = {};
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        cv(end+1) = std(relFP(strcmp(valid_rxns,rxnLabel{j}),:))/mean(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end
            
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & cv > 0.1)>0));


figure;
myLevel = relFP(strcmp(valid_rxns,rxnID1),:);
fit = fitlm(myLevel,myFlux);
h = plot(fit);
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','o')
set(h(1),'MarkerSize',7)
set(h(2), {'color'},{'k'}) 
set(h(3), {'color'},{'k'}) 
set(h(4), {'color'},{'k'}) 
xlabel('Relative Flux Potential');
ylabel('Relative Flux (absolute value)');
text(0.6,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
text(0.6,max(myFlux)*0.85,['FDR = ',num2str(fdr_r(strcmp(testedRxn,rxnID1)),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
% plt.XTick = -1:0.2:1;
% plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Title = rxnID1;
plt.Interpreter = 'None';
plt.export(['figures/others/new_rFP_expression_',rxnID1,'.tiff']);

%% check expression consistency - discountinued
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
