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

fdr_rho = mafdr(p_rho,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))

fprintf('%d rxns give significant positive correlation by spearman\n',sum(rho(fdr_rho<0.05)>0));
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

%% save the data for R plot
Candidates = testedRxn;
candidatesP = r;
[~,A] = sort(candidatesP,'descend');
Candidates = Candidates(A);
flux_out = [];
exp_out = [];
%j = 15;
for j = 1:length(Candidates)  
    fluxMeasure = abs(fluxMat(strcmp(rxnLabel,Candidates{j}),:));
    expression = normalizedLevel(strcmp(valid_rxns,Candidates{j}),:);
    flux_out = [flux_out, fluxMeasure'];
    exp_out = [exp_out, expression'];
end
t = array2table(flux_out);
t.Properties.VariableNames = Candidates;
writetable(t,'output/flux_expression_correlation_flux_mat_to_plot.csv','WriteRowNames',1);
t = array2table(exp_out);
t.Properties.VariableNames = Candidates;
writetable(t,'output/flux_expression_correlation_exp_mat_to_plot.csv','WriteRowNames',1);



