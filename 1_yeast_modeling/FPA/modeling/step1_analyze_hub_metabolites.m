%% About
% analyze the correlation between flux and correlation, and generates the
% rxn-centered expression levels for all downstream analyses 
%% load the model
addpath('./../scripts/')
model = loadYeatModel();
%% load the expression files, distance matrix, and other optional inputs
% load proteomics (the original raw proteomics data in SIMMER paper)
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
% we use the non-log level
proTbl{:,2:end} = 2.^proTbl{:,2:end};
% load the biomass composition in SIMMER paper
coeffTbl = readtable('./../input/YeastJoshua/originalDataTbl/coefficients_gram_per_gDW.xlsx');
% preprocess the expression table
% the FPA matrix will be in the same order as the master_expression
conditions = proTbl.Properties.VariableNames(2:end);
master_expression_perPro = {};
master_expression_perDW= {};
geneInd = ismember(proTbl.Gene, model.genes); % get the index of genes in the model
for i = 1:length(conditions)
    expression = struct();
    expression.genes = proTbl.Gene(geneInd);
    expression.value = proTbl.(conditions{i})(geneInd);
    master_expression_perPro{i} = expression;
    expression.value = expression.value * coeffTbl.(conditions{i})(strcmp(coeffTbl.Metabolite,'Protein'));
    master_expression_perDW{i} = expression;
end
master_expression_pro_perPro = master_expression_perPro;% this is the reletive level
master_expression_pro_perDW = master_expression_perDW;% this is the absolute level
load('output/normalizedLevels_partialExcluded.mat')
%% correlation with ROI expression: relative levels (unit free)
normalizedLevel = normalizedLevel_pro_perPro;
valid_rxns = valid_rxns_pro_perPro ;
% flux table 
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
fluxMat = [];
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
% normalize flux unit
% when normalize the flux by the flux of biomass production (growth rate),
% the unit of growth rate needs to be taken care of. In chemostat setting,
% steady state was defined as stable OD (see SIMMER paper), which means
% steady cell density (number, aka, volume). Therefore, the dilution rate
% is a measure of per cell flux. So, we should normalize the internal flux
% under /ml cell metric
% flux is in  (moles / hr / mL cells); no conversion is needed. 

GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat = fluxMat ./ repmat(GRrate.DR_Actual',size(fluxMat,1),1);
fluxMat_normalized = fluxMat;
%Computing the correlation
rho=[];
p_rho=[];
r=[];
p_r=[];
testedRxn = {};
rxnLabel = fluxTbl.Model_Reaction_ID;
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        [rho(end+1),p_rho(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
        [r(end+1),p_r(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
    end
end
fdr_rho = mafdr(p_rho,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))

fprintf('%d rxns give significant positive correlation by spearman\n',sum(rho(fdr_rho<0.05)>0));
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

%% analyze the hub metabolite relation
% hubness = [];
% medianCorr = [];
% porpSigCorr = [];
% for i = 1:length(model.mets)
%     hubness(i) = sum(model.S(i,:)~=0);% total number of connections (indegree + outdegree)
%     assoRxns = model.rxns(model.S(i,:)~=0);
%     if isempty(intersect(assoRxns, testedRxn))
%         medianCorr(i) = NaN;
%         porpSigCorr(i) = NaN;
%     else
%         corrVals = r(ismember(testedRxn,assoRxns));
%         fdrVals = fdr_r(ismember(testedRxn,assoRxns));
%         medianCorr(i) = median(corrVals);
%         porpSigCorr(i) = sum(fdrVals < 0.05 & corrVals > 0) / length(corrVals);
%     end
% end
% plot(log2(hubness), medianCorr,'.');
% plot(log2(hubness), porpSigCorr,'.');

% we check the mutual informing rate for all the two-reaction unit
% and its relation to the hubness of the bridging metabolite 

% load distance matrix to save time on calculating the hubness
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

% make the hubness matrix 
% first get all two-reaction modules from the orginal distance matrix
cMat_tent = full(distMat_normal <= 1);
labels_ID = regexprep(labels,'_.$','');
a = zeros(length(model.rxns),length(model.rxns));
inds = false(length(model.rxns),length(labels_ID));
for j = 1:length(model.rxns)
    inds(j,:) = strcmp(labels_ID, model.rxns{j});
end
for i = 1:length(model.rxns)
    ind = inds(i,:);
    if any(ind)
        tmp = max(cMat_tent(ind,:),[],1);
        for j = 1:length(model.rxns)
            ind = inds(j,:);
            if any(ind)
                a(i,j) = max(tmp(ind));
            end
        end
    end
    i
end
cMat = a; % 0 or NaN: not connected, 1: connected 

% next get the hubness (wtd distance * 4) from the wtd distance 
hubMat_tent = distMat_wtd;
labels_ID = regexprep(labels,'_.$','');
a = zeros(length(model.rxns),length(model.rxns));
inds = false(length(model.rxns),length(labels_ID));
for j = 1:length(model.rxns)
    inds(j,:) = strcmp(labels_ID, model.rxns{j});
end
for i = 1:length(model.rxns)
    ind = inds(i,:);
    if any(ind)
        tmp = min(hubMat_tent(ind,:),[],1);
        for j = 1:length(model.rxns)
            ind = inds(j,:);
            if any(ind)
                a(i,j) = min(tmp(ind));
            end
        end
    end
    i
end
hubnessMat = cMat .* a .* 4; % number of connections
t = array2table(hubnessMat);
t.Properties.VariableNames = model.rxns;
t.Properties.RowNames = model.rxns;
writetable(t,'output/NoTrack_hubnessMat.csv','WriteRowNames',1);

% next prepare the mutually informing matrix (whether rxnA's expression
% dictates the flux of rxnB)
mutInfoMat = NaN(length(model.rxns), length(model.rxns));
ExpMeasured = ismember(model.rxns,valid_rxns);
FluxMeasured = ismember(model.rxns,rxnLabel);
for i = 1:length(model.rxns)
    for j = 1:length(model.rxns)
        mutinformed = [];
        if FluxMeasured(i) && ExpMeasured(j)
            % calculate
            fluxMeasure = fluxMat(strcmp(rxnLabel,model.rxns{i}),:);
            expressionMeasure = normalizedLevel(strcmp(valid_rxns,model.rxns{j}),:);
            [r0,p_r0] = corr(expressionMeasure',abs(fluxMeasure)','type','Pearson');
            mutinformed = [mutinformed,-log10(p_r0) * sign(r0)];% we save both pvalue and the direction of correlation
        end
        
        if FluxMeasured(j) && ExpMeasured(i)
            % calculate
            fluxMeasure = fluxMat(strcmp(rxnLabel,model.rxns{j}),:);
            expressionMeasure = normalizedLevel(strcmp(valid_rxns,model.rxns{i}),:);
            [r0,p_r0] = corr(expressionMeasure',abs(fluxMeasure)','type','Pearson');
            mutinformed = [mutinformed,-log10(p_r0) * sign(r0)];
        end
        if ~isempty(mutinformed)
            mutInfoMat(i,j) = max(mutinformed);% the most significant positive correlation we got
        end
    end
    i
end
t = array2table(mutInfoMat);
t.Properties.VariableNames = model.rxns;
t.Properties.RowNames = model.rxns;
writetable(t,'output/NoTrack_mutInfoMat.csv','WriteRowNames',1);

% compare the hubness vs. mutinfo probalilty for all two-reaction modules
rxnA = {};
rxnB = {};
hubness = [];
mutInformed = [];
for i = 1:(length(model.rxns)-1)
    for j = (i+1):length(model.rxns) % all unique two-reaction modules
        if hubnessMat(i,j) > 1 && ~isnan(mutInfoMat(i,j))% connected and measured
            rxnA = [rxnA; model.rxns(i)];
            rxnB = [rxnB; model.rxns(j)];
            hubness = [hubness; hubnessMat(i,j)];
            mutInformed = [mutInformed; mutInfoMat(i,j)];
        end
    end
    i
end
% calculate for fdr 
mutInformed = sign(mutInformed) .* mafdr(10.^(-abs(mutInformed)),'BHFDR', true);
% we consider any reaction pair with significant positive correlation (fdr
% < 0.05) as mutually informed
mutInformed = mutInformed > 0 & mutInformed < 0.05; 
%% calculate and plot 
ValidHubnessVals = sort(unique(hubness)); 
% counts = tabulate(hubness);
informRate = []; % for each hubness of bridging metabolite, how many percents of two-reaction modules were mutually informed
dataCount = [];
for i = 1:length(ValidHubnessVals)
    informRate(i) = sum(mutInformed(hubness == ValidHubnessVals(i))) / sum(hubness == ValidHubnessVals(i));
    dataCount(i) = sum(hubness == ValidHubnessVals(i));
end
% we group the hubness when the sample size is too small (since to total
% measurement are only 232 fluxes and ~600 expressions. not all
% two-reaction modules were quantified)
minSample = 20;
clean = false;
while ~clean
    isend = false;
    i = 1;
    while ~isend
        if dataCount(i) < minSample
            % group with the next hubness val
            % the new hubness represent the average hubness and new
            % infoRate is the avergae after we group all these two-reaction
            % modules
            informRate(i+1) = (informRate(i) .* dataCount(i) + informRate(i+1) .* dataCount(i+1))/(dataCount(i) + dataCount(i+1));
            ValidHubnessVals(i+1) = (ValidHubnessVals(i) .* dataCount(i) + ValidHubnessVals(i+1) .* dataCount(i+1))/(dataCount(i) + dataCount(i+1));
            dataCount(i+1) = dataCount(i) + dataCount(i+1);
            informRate(i) = [];
            ValidHubnessVals(i) = [];
            dataCount(i) = [];
        else
            i = i+1;
        end
        isend = i == length(dataCount);
    end
    clean = ~any(dataCount < minSample);
end
median(dataCount)
% plot(log2(ValidHubnessVals), informRate,'.');
informRate = informRate .* 100; % convert to percentages for plotting
%% save plot
fit = polyfit(log2(ValidHubnessVals),informRate,4);

x1 = linspace(0.9,0.1+max(log2(ValidHubnessVals)));
y1 = polyval(fit,x1);
figure
h1 = plot(log2(ValidHubnessVals), informRate,'o');
hold on
h2 = plot(x1,y1);
hold off
set(h1, {'color'},{'k'}) 
set(h1,'MarkerSize',5)
set(h1,'LineWidth',1)
set(h2,'LineWidth',1)
xlabel('Degree of the bridging metabolite (log2)');
ylabel('Cross-informing rate (%)');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.export(['figures/two_reaction_module_mutually_informing_rate.pdf']);

