%% rationale
% since the flux clustered into coflux modules. one questions is whether
% the expression of genes in purine module could correlate with less
% related module (argine and histidine). If so, it indicates that the
% transcription pattern is something shared by the big cluster; if not, it
% inidicates individual pathway is fine-tuned. 

% ==> see figure in R; we found that puring and lysine is a module; while
% other AA is a module (his, phe, his)
%% load the model
addpath('./../scripts/')
model = loadYeatModel();
%% load the expression files, distance matrix, and other optional inputs
% load proteomics
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
% we use the non-log level
proTbl{:,2:end} = 2.^proTbl{:,2:end};
% normalizaed
coeffTbl = readtable('./../input/YeastJoshua/originalDataTbl/coefficients_gram_per_gDW.xlsx');
% preprocess the expression table
% to facilate the future use of the expression of many samples, we
% re-organize it into a structure variable.
% the FPA matrix will be in the same order as the master_expression
conditions = proTbl.Properties.VariableNames(2:end);
% make a new master_expression for these four conditions.
master_expression_perPro = {};% we call this variable "master_expression"
master_expression_perDW= {};% we call this variable "master_expression"
geneInd = ismember(proTbl.Gene, model.genes); % get the index of genes in the model
for i = 1:length(conditions)
    expression = struct();
    expression.genes = proTbl.Gene(geneInd);
    expression.value = proTbl.(conditions{i})(geneInd);
    master_expression_perPro{i} = expression;
    expression.value = expression.value * coeffTbl.(conditions{i})(strcmp(coeffTbl.Metabolite,'Protein'));
    master_expression_perDW{i} = expression;
end
master_expression_pro_perPro = master_expression_perPro;
master_expression_pro_perDW = master_expression_perDW;
%%
load('output/normalizedLevels_partialExcluded.mat')
%% correlation with immediate vincinity: protein : normalized unit free
close all
normalizedLevel = normalizedLevel_pro_perPro;
valid_rxns = valid_rxns_pro_perPro ;
% load pathway info
pathways = readtable('manualPathwayAnnotation.csv');
% purine expression 
purineExp = normalizedLevel(ismember(valid_rxns,pathways.rxn(strcmp(pathways.pathway_major,'Purine metabolism'))),:);
purineRxns = valid_rxns(ismember(valid_rxns,pathways.rxn(strcmp(pathways.pathway_major,'Purine metabolism'))));

% flux table 
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
fluxMat = [];
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat = fluxMat ./ repmat(GRrate.DR_Actual',size(fluxMat,1),1);
fluxMat_normalized = fluxMat;
rxnLabel = fluxTbl.Model_Reaction_ID;
% argine and histidine
argFlux = fluxMat_normalized(ismember(rxnLabel,pathways.rxn(strcmp(pathways.pathway_major,'Arginine biosynthesis'))),:);
argRxns = rxnLabel(ismember(rxnLabel,pathways.rxn(strcmp(pathways.pathway_major,'Arginine biosynthesis'))));


% control
r0=zeros(length(purineRxns),length(purineRxns));
p_r0=ones(length(purineRxns),length(purineRxns));
for j = 1:length(purineRxns)
    myExp = normalizedLevel(strcmp(valid_rxns,purineRxns{j}),:)';
    for i = 1:length(purineRxns)
        myFlux = abs(fluxMat(strcmp(rxnLabel,purineRxns{i}),:))';
        [a1,a2] = corr(myExp,myFlux,'type','Pearson');
        r0(j,i) = a1;
        p_r0(j,i) = a2;
    end
end


% compute correlation: purine expression vs. arg flux
r=zeros(length(purineRxns),length(argRxns));
p_r=ones(length(purineRxns),length(argRxns));
for j = 1:length(purineRxns)
    myExp = normalizedLevel(strcmp(valid_rxns,purineRxns{j}),:)';
    for i = 1:length(argRxns)
        myFlux = abs(fluxMat(strcmp(rxnLabel,argRxns{i}),:))';
        [a1,a2] = corr(myExp,myFlux,'type','Pearson');
        r(j,i) = a1;
        p_r(j,i) = a2;
    end
end

HeatMap(-log10([p_r0,p_r]))

%% Computing the correlation between any expression/reaction pair 
myrxns = intersect(valid_rxns,rxnLabel);

r=zeros(length(myrxns),length(myrxns));
p_r=ones(length(myrxns),length(myrxns));
for j = 1:length(myrxns)
    myExp = normalizedLevel(strcmp(valid_rxns,myrxns{j}),:)';
    for i = 1:length(myrxns)
        myFlux = abs(fluxMat(strcmp(rxnLabel,myrxns{i}),:))';
        [a1,a2] = corr(myExp,myFlux,'type','Pearson');
        r(j,i) = a1;
        p_r(j,i) = a2;
    end
end

clustergram(-log10([p_r]))
% save data 
t = array2table(r);
t.Properties.RowNames = myrxns;
writetable(t,'output/crossCorrelation_rMat.csv','WriteRowNames',1);
t = array2table(p_r);
t.Properties.RowNames = myrxns;
writetable(t,'output/crossCorrelation_pMat.csv','WriteRowNames',1);






