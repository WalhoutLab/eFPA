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
%% make the rxn-centered expression levels
parpool(4)
%% calculate penalty (normalized level) - reletive expression
penalty = calculatePenalty_partialExcluded(model,master_expression_pro_perPro);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only valid (has data) level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns = model.rxns(~rmInd);
normalizedLevel_pro_perPro = normalizedLevel;
valid_rxns_pro_perPro = valid_rxns;
%% calculate penalty (normalized level) - absolute expression
penalty = calculatePenalty_partialExcluded(model,master_expression_pro_perDW);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only valid (has data) level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns = model.rxns(~rmInd);
normalizedLevel_pro_perDW = normalizedLevel;
valid_rxns_pro_perDW = valid_rxns;
save('output/normalizedLevels_partialExcluded.mat','normalizedLevel_pro_perPro','valid_rxns_pro_perPro','normalizedLevel_pro_perDW','valid_rxns_pro_perDW');
%%
load('output/normalizedLevels_partialExcluded.mat')
%% correlation with ROI expression: absolute measures
normalizedLevel = normalizedLevel_pro_perDW;
valid_rxns = valid_rxns_pro_perDW ;
% flux table (raw flux in SIMMER dataset. the rxnID was updated based on
% the used yeast model)
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% use the mid point of FVA interval as the flux values
fluxMat =[];
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
% normalize flux unit
% flux is in  (moles / hr / mL cells); could be further normalized to
% mmole/hr/gDW by chemostat info: gDCW/ml. such that it is comparable with
% the per gDW protein content
dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
factor = repmat(dwTbl.gDCW_mL',size(fluxMat,1),1);
fluxMat = fluxMat * 1000 ./ factor; %mmoles/hr/gDW
fluxMat_raw = fluxMat;

%Computing the correlation between expression and flux
rho=[];
r = [];
p_rho=[];
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
% correct multiple testing 
fdr_rho = mafdr(p_rho,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))

fprintf('%d rxns give significant positive correlation by spearman\n',sum(rho(fdr_rho<0.05)>0));
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

sigCorrRxns = testedRxn(fdr_r<0.05 & r >0);
corrTbl = table(testedRxn',r',fdr_r');
save('output/sig_corr_rxns_raw2raw.mat','sigCorrRxns','corrTbl');

figure(1)
hold on
histogram(rho,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(rho(fdr_rho<0.05 & rho >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
legend({'all testable reactions','significantly correlated reactions'})
hold off
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export('figures/Expression_flux_correlation_rawFlux_spearman.pdf');

figure(2)
hold on
histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(r(fdr_r<0.05 & r >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
hold off
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export('figures/Expression_flux_correlation_rawFlux_pearson.pdf');
%% correlation with ROI expression: relative levels (unit free)
close all
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
% save result
t = array2table(r');
t.Properties.RowNames = testedRxn;
writetable(t,'output/flux_expression_pearson_correlation_rel2rel.csv','WriteRowNames',1);

fdr_rho = mafdr(p_rho,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))

fprintf('%d rxns give significant positive correlation by spearman\n',sum(rho(fdr_rho<0.05)>0));
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

sigCorrRxns = testedRxn(fdr_r<0.05 & r >0);
corrTbl = table(testedRxn',r',fdr_r');
save('output/sig_corr_rxns_rel2rel.mat','sigCorrRxns','corrTbl');

% annotate the significant reactions
[ROI_annotated, ROI_enrichment] = annotateRxnSet(sigCorrRxns,model,testedRxn);
ROI_annotated(:,5) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl.Var1);
ROI_annotated(A,5) = mat2cell(corrTbl{B(A),2},ones(sum(A),1),1); % pearson r
ROI_annotated(:,6) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl.Var1);
ROI_annotated(A,6) = mat2cell(corrTbl{B(A),3},ones(sum(A),1),1); % FDR
ROI_annotated = [{'rxnID','formula','GPR','firstSubsystem','pearson r','FDR'};ROI_annotated];
ROI_enrichment = [{'SubSystems','number of significantly correlated rxns','number of total rxns in the subsys','FDR'};ROI_enrichment];
ROI_annotated_pos = ROI_annotated;
ROI_enrichment_pos = ROI_enrichment;

% annotate the significant reactions - the negative correlations
sigCorrRxns_neg = testedRxn(fdr_r<0.05 & r < 0);
[ROI_annotated, ROI_enrichment] = annotateRxnSet(sigCorrRxns_neg,model,testedRxn);
ROI_annotated(:,5) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl.Var1);
ROI_annotated(A,5) = mat2cell(corrTbl{B(A),2},ones(sum(A),1),1); % pearson r
ROI_annotated(:,6) = repmat({''},size(ROI_annotated,1),1);
[A B] = ismember(ROI_annotated(:,1),corrTbl.Var1);
ROI_annotated(A,6) = mat2cell(corrTbl{B(A),3},ones(sum(A),1),1); % FDR
ROI_annotated = [{'rxnID','formula','GPR','firstSubsystem','pearson r','FDR'};ROI_annotated];
ROI_enrichment = [{'SubSystems','number of significantly correlated rxns','number of total rxns in the subsys','FDR'};ROI_enrichment];
ROI_annotated_neg = ROI_annotated;
ROI_enrichment_neg = ROI_enrichment;

figure(4)
hold on
histogram(rho,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(rho(fdr_rho<0.05 & rho >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
hold off
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export('figures/Expression_flux_correlation_relFlux_spearman.pdf');

figure(5)
hold on
histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(r(fdr_r<0.05 & r >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
hold off
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export('figures/Expression_flux_correlation_relFlux_pearson.pdf');
%% a special look at the central carbon metabolism by their correlation
PathwayAnn = readtable('./pathway_annotations.xlsx','Sheet','more_info');
centralCarbon = {'Dihydroorotate production','Glycolysis','Methylglyoxal metabolism','Oxidative phosphorylation',...
    'Pentose phosphate pathway','Pyruvate metabolism','TCA cycle'};
for i = 1:length(centralCarbon)
    myRxnSet = PathwayAnn.rxn(strcmp(PathwayAnn.manual_pathway,centralCarbon{i}));
    figure;
    hold on
    histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
    xlim([-1,1]);
    xlabel('Correlation coefficient');
    ylabel('Number of reactions');
    histogram(r(ismember(testedRxn, myRxnSet)),'FaceColor','#D95319','BinEdges',-1:0.2:1)
    legend({'all testable reactions',sprintf('%s',centralCarbon{i})})
    hold off
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [1.95, 1.6125];
    plt.LineWidth = 1;
    plt.FontSize = 10;
    plt.XTick = -1:0.2:1;
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.ShowBox = 'off';
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.export(['figures/centralCarbon_example/Expression_flux_correlation_relFlux_pearson_',centralCarbon{i},'.pdf']);
end
myRxnSet = PathwayAnn.rxn(ismember(PathwayAnn.manual_pathway,centralCarbon));
figure;
hold on
histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(r(ismember(testedRxn, myRxnSet)),'FaceColor','#D95319','BinEdges',-1:0.2:1)
legend({'all testable reactions',sprintf('central carbon reactions')})
hold off
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/centralCarbon_example/Expression_flux_correlation_relFlux_pearson_centralCarbon.pdf']);
%% individual scattering plot -- save all plots 
close all
Candidates = testedRxn;
%j = 15;
for j = 1:length(Candidates)  
    fluxMeasure = fluxMat(strcmp(rxnLabel,Candidates{j}),:);
    expression = normalizedLevel(strcmp(valid_rxns,Candidates{j}),:);
    
    fluxMeasure_raw = fluxMat_raw(strcmp(rxnLabel,Candidates{j}),:);
    expression_raw = normalizedLevel_pro_perDW(strcmp(valid_rxns_pro_perDW,Candidates{j}),:);

    fit_rel = fitlm(expression,abs(fluxMeasure));
    
    figure(4)
    h = plot(fit_rel);
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
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [2, 1.75];
    plt.LineWidth = 1;
    plt.FontSize = 10;
    plt.XTick = -1:0.2:1;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    title([Candidates{j},' r = ',num2str(sqrt(fit_rel.Rsquared.Ordinary))],'Interpreter','none');
    plt.ShowBox = 'off';
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.TickDir = 'out';
    plt.export(['figures/correlationPlots/',Candidates{j},'_rel2rel.pdf']);
    j
end
%% enrichment of direct connections
byProducts = {'carbon dioxide';'AMP';'NADP(+)';'NADPH';'diphosphate';'oxygen';'NADH';'NAD';'phosphate';'ADP';'coenzyme A';'ATP';'H2O';'H+';'GTP';'GDP';'(R)-carnitine';'FAD';'FADH2'};
% add compartment label to byproducts
byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,'\s\[.*\]$',''),model.metNames, 'UniformOutput',false),byProducts));
model_bak = model;
model.S(ismember(model.mets,byProducts),:) = [];
model.mets(ismember(model.mets,byProducts)) = [];

cMat = zeros(length(model.rxns),length(model.rxns));% connection matrix
for i = 1:length(model.rxns)
    % out degree
    if model.lb(i) < -1e-10
        metsID = model.S(:,i)<0;
        nextRxnID = any(model.S(metsID,:),1);
        nextRxnID(i) = false;
        rxnset = model.rxns(nextRxnID);
        r1 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)>0,1)');...
            rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)<0,1)')];
    else
        r1 = {};
    end
    if model.ub(i) > 1e-10
        metsID = model.S(:,i)>0;
        nextRxnID = any(model.S(metsID,:),1);
        nextRxnID(i) = false;
        rxnset = model.rxns(nextRxnID);
        r2 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)>0,1)');...
            rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)<0,1)')];
    else
        r2 = {};
    end
    
    % in degree
    if model.lb(i) < -1e-10
        metsID = model.S(:,i)>0;
        nextRxnID = any(model.S(metsID,:),1);
        nextRxnID(i) = false;
        rxnset = model.rxns(nextRxnID);
        r3 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)<0,1)');...
            rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)>0,1)')];
    else
        r3 = {};
    end
    if model.ub(i) > 1e-10
        metsID = model.S(:,i)<0;
        nextRxnID = any(model.S(metsID,:),1);
        nextRxnID(i) = false;
        rxnset = model.rxns(nextRxnID);
        r4 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)<0,1)');...
            rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)>0,1)')];
    else
        r4 = {};
    end
    r_c = unique([r1;r2;r3;r4]);
    cMat(i,ismember(model.rxns,r_c)) = 1;
end
model = model_bak;       

% analysis 
Candidates = testedRxn(fdr_r<0.05 & r >0);
candidatesP = fdr_r(fdr_r<0.05 & r >0);
[~,A] = sort(candidatesP);
Candidates = Candidates(A);
hitRxns = Candidates;
allRxnSet_tested = testedRxn;
allRxnSet_GPR = model.rxns(cellfun(@(x) ~isempty(x),model.grRules));

c_hit = sum(sum(cMat(ismember(model.rxns,hitRxns),ismember(model.rxns,hitRxns))));
% permutation test - against all tested rxns
c_perm = zeros(100000,1);
rng(1126)
for i = 1: 100000
    sample = allRxnSet_tested(randperm(length(allRxnSet_tested),length(hitRxns)));
    c_perm(i) = sum(sum(cMat(ismember(model.rxns,sample),ismember(model.rxns,sample))));
end

figure;
histogram(c_perm,'Normalization','pdf','BinEdges',0:2:140)
xline(c_hit,'r')
xlabel('Total number of connections','FontSize', 12)
ylabel('Frequency','FontSize', 12)
txt = ['p = ', num2str((1+ sum(c_perm>=c_hit)) / (1+length(c_perm)))];
text(c_hit+5,0.025,txt,'HorizontalAlignment','left','FontSize', 10)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.25];
plt.LineWidth = 0.75;
plt.FontSize = 10;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.export('figures/Expression_flux_correlation_connectivityEnrichment.pdf');
