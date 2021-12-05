%% About
% calculate and save the coexpression (PCC) tables
%% load inputs
% condition names
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);
% flux tables
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
% load expression 
load('output/normalizedLevels_partialExcluded.mat');
%% correlation of all expressions 
myrxns = intersect(valid_rxns,rxnLabel);
normalizedLevel = normalizedLevel_pro_perPro;
valid_rxns = valid_rxns_pro_perPro ;
myExp = normalizedLevel(ismember(valid_rxns,myrxns),:);
myExpID = valid_rxns(ismember(valid_rxns,myrxns));
[cor, p] = corr(myExp');
clustergram(cor,'RowLabels',myExpID)
t = array2table(myExp);
t.Properties.RowNames = myExpID;
t.Properties.VariableNames = conditions;
writetable(t,'output/relativeExpressionTable.csv','WriteRowNames',1);

%% annotate 
ROI = mygroup1.RowNodeNames;
[ROI_annotated, ROI_enrichment] = annotateRxnSet(ROI,model,testedRxn);

%% search TF for rxn clusters and see what is measured
TFmat = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202128_expressionAndBinding.txt',...
    'Delimiter',';','HeaderLines',1);
TargetgeneName = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202128_expressionAndBinding.txt',...
    'Delimiter',';','ReadVariableNames',0);
targetGeneName = TargetgeneName{1,:};
TFnames = TFmat.Var1;
TFmat = TFmat{:,ismember(targetGeneName,model.geneNames)};
targetGeneName = targetGeneName(ismember(targetGeneName,model.geneNames));
% CBF1
% cbf1 expression
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
cbf1 = proTbl(strcmp(proTbl.Gene,'YJR060W'),:);
cbf1{:,2:end} = 2.^cbf1{:,2:end};
writetable(cbf1,'output/CBF1_expression.csv','WriteRowNames',1);
% cbf1 regulated rxns
cbf1Rxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Cbf1p'),:)==1))),2)));
writetable(cbf1Rxns,'output/CBF1_targetRxn.csv','WriteRowNames',1);

%% for the few measured TF, see if some is correlated 
TFmat = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202128_expressionAndBinding_measuredOnly.txt',...
    'Delimiter',';','HeaderLines',1);
TargetgeneName = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202128_expressionAndBinding_measuredOnly.txt',...
    'Delimiter',';','ReadVariableNames',0);
targetGeneName = TargetgeneName{1,:};
TFnames = TFmat.Var1;
TFmat = TFmat{:,ismember(targetGeneName,model.geneNames)};
% the measured TF
TFnames(sum(TFmat,2)>0)
% regulated rxns
cbf1Rxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Cbf1p'),:)==1))),2)));
Tup1Rxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Tup1p'),:)==1))),2)));
Rap1Rxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Rap1p'),:)==1))),2)));

targetGeneName = targetGeneName(ismember(targetGeneName,model.geneNames));
% Tup1 expression
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
Tup1 = proTbl(strcmp(proTbl.Gene,'YCR084C'),:);
Tup1{:,2:end} = 2.^Tup1{:,2:end};
writetable(Tup1,'output/Tup1_expression.csv','WriteRowNames',1);
% Tup1 regulated rxns
Tup1Rxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Tup1p'),:)==1))),2)));
writetable(Tup1Rxns,'output/Tup1_targetRxn.csv','WriteRowNames',1);

% Rap1 expression
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
Rap1 = proTbl(strcmp(proTbl.Gene,'YNL216W'),:);
Rap1{:,2:end} = 2.^Rap1{:,2:end};
writetable(Rap1,'output/Rap1_expression.csv','WriteRowNames',1);
% Rap1 regulated rxns
Rap1Rxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Rap1p'),:)==1))),2)));
writetable(Rap1Rxns,'output/Rap1_targetRxn.csv','WriteRowNames',1);


% add a TF once a time to see a overall correlation with a cluster 




%% check if some TF covers a precise cluster 
[cor, p] = corr(fluxMat_normalized');
clustergram(cor,'RowLabels',targetRxns);

geneset = model.genes(any(model.rxnGeneMat(ismember(model.rxns,mygroup.RowNodeNames),:)));

TFmat = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202129_expressionPlusBinding.txt',...
    'Delimiter',';','HeaderLines',1);
TargetgeneName = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202129_expressionPlusBinding.txt',...
    'Delimiter',';','ReadVariableNames',0);
targetGeneName = TargetgeneName{1,:};
TFnames = TFmat.Var1;
TFmat = TFmat{:,ismember(targetGeneName,model.geneNames)};
targetGeneName = targetGeneName(ismember(targetGeneName,model.geneNames));


% calculate the enrichment of specifically regulating a set of target gene 
p_enri = [];
allGene = model.geneNames;%(any(model.rxnGeneMat(ismember(model.rxns,targetRxns),:)));
geneset2 = model.geneNames(ismember(model.genes,geneset));
for i = 1:size(TFmat,1)
    reg = targetGeneName(TFmat(i,:) == 1);
    N_in = intersect(reg, geneset2);
    N_all = intersect(reg, allGene);
    p_enri(i,1) = 1- hygecdf(length(N_in)-1,length(allGene),length(N_all),length(geneset));
end
TFnames_sorted = TFnames;
[~,I] = sort(p_enri);
TFnames_sorted = TFnames_sorted(I);
TFnames_sorted(1:10)


% Gcn4p
% Gcn4p regulated rxns
Gcn4pRxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Gcn4p'),:)==1))),2)));
writetable(Gcn4pRxns,'output/Gcn4_targetRxn.csv','WriteRowNames',1);

% Lys14p
Lys14pRxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Lys14p'),:)==1))),2)));
writetable(Lys14pRxns,'output/Lys14_targetRxn.csv','WriteRowNames',1);

% Leu3p
Leu3pRxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Leu3p'),:)==1))),2)));
writetable(Leu3pRxns,'output/Leu3_targetRxn.csv','WriteRowNames',1);

% Gln3p
Gln3pRxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Gln3p'),:)==1))),2)));
writetable(Gln3pRxns,'output/Gln3_targetRxn.csv','WriteRowNames',1);

% Gcr2p?/ ==> no specific regulator found; all checked are master regualtor
% with enrichment
Gcr2pRxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName(TFmat(strcmp(TFnames,'Rap1p'),:)==1))),2)));
writetable(Gln3pRxns,'output/Gln3_targetRxn.csv','WriteRowNames',1);


%%
TFmat(all(TFmat==0,2),:)= [];
dist= squareform(pdist(TFmat','cosine'));
% dist= corr(TFmat);
clustergram(dist,'RowLabels',targetGeneName);


%% co-expression of TF and metabolic reactions based on microarray data 
%% load the model
addpath('./../scripts/')
model = loadYeatModel();
%% load the expression files, distance matrix, and other optional inputs
% load microarray ==> we dont analyze RNA data since it is not clean
rnaTbl = readtable('./../input/YeastJoshua/MicroArray/matched_knn_imputed_log2_FC_to_reference_pmid_17959824.txt');% this is the log2(FC_reference)
% we use the non-log level
rnaTbl{:,3:end} = 2.^rnaTbl{:,3:end};
% save
writetable(rnaTbl,'output/RNA_all_expression_FC_table.csv','WriteRowNames',0,'QuoteStrings',1);
coeffTbl = readtable('./../input/YeastJoshua/originalDataTbl/coefficients_gram_per_gDW.xlsx');
conditions = rnaTbl.Properties.VariableNames(3:end);
% preprocess the expression table
% to facilate the future use of the expression of many samples, we
% re-organize it into a structure variable.
% the FPA matrix will be in the same order as the master_expression
% make a new master_expression for these four conditions.
master_expression_perRNA = {};% we call this variable "master_expression"
master_expression_perDW= {};% we call this variable "master_expression"
geneInd = ismember(rnaTbl.YORF, model.genes); % get the index of genes in the model
for i = 1:length(conditions)
    expression = struct();
    expression.genes = rnaTbl.YORF(geneInd);
    expression.value = rnaTbl.(conditions{i})(geneInd);
    master_expression_perRNA{i} = expression;
    expression.value = expression.value * coeffTbl.(conditions{i})(strcmp(coeffTbl.Metabolite,'RNA'));
    master_expression_perDW{i} = expression;
end
master_expression_rna = master_expression_perRNA;
master_expression_rna_perDW = master_expression_perDW;
%% parpool
parpool(4)
%% calculate penalty (normalized level)
penalty = calculatePenalty(model,master_expression_rna);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only data-derived level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns = model.rxns(~rmInd);
% save the normalized level table
RNAtable = array2table(normalizedLevel);
RNAtable.Properties.RowNames = valid_rxns;
RNAtable.Properties.VariableNames = conditions;
writetable(RNAtable,'output/RNA_rxn_expression_FC_table.csv','WriteRowNames',1);


%% save the regulation matrix related to metabolism 
TFmat = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202129_expressionPlusBinding.txt',...
    'Delimiter',';','HeaderLines',1);
TargetgeneName = readtable('./../input/YeastJoshua/RegulationMatrix_Documented_202129_expressionPlusBinding.txt',...
    'Delimiter',';','ReadVariableNames',0);
targetGeneName = TargetgeneName{1,:};
TFnames = TFmat.Var1;
TFmat = TFmat{:,ismember(targetGeneName,model.geneNames)};
targetGeneName = targetGeneName(ismember(targetGeneName,model.geneNames));

rmInd = all(TFmat==0,2);
TFmat(rmInd,:) = [];
TFnames(rmInd) = [];

% make the rxn level GRN
rxnGRN = zeros(length(TFnames),length(model.rxns));
for i = 1:length(model.rxns)
    if any(ismember(targetGeneName,model.geneNames(logical(model.rxnGeneMat(i,:)))))
        rxnGRN(:,i) = max(TFmat(:,ismember(targetGeneName,model.geneNames(logical(model.rxnGeneMat(i,:))))),[],2);
    end
end

TFmat= array2table(TFmat);
rxnGRN = array2table(rxnGRN);
TFnames = table(TFnames);
targetRxnNames = table(model.rxns);
targetGeneName = table(targetGeneName');
writetable(TFnames,'output/regulationMatrix_metabolicOnly_TFnames.csv','WriteRowNames',0);
writetable(targetGeneName,'output/regulationMatrix_metabolicOnly_targetGeneName.csv','WriteRowNames',0);
writetable(TFmat,'output/regulationMatrix_metabolicOnly_matrix.csv','WriteRowNames',0);
writetable(rxnGRN,'output/regulationMatrix_metabolicOnly_matrix_rxn.csv','WriteRowNames',0);
writetable(targetRxnNames,'output/regulationMatrix_metabolicOnly_targetRxnNames.csv','WriteRowNames',0);

%% TARGET genes
YHP1pRxns = table(model.rxns(any(model.rxnGeneMat(:,ismember(model.geneNames,targetGeneName.Var1(TFmat{strcmp(TFnames.TFnames,'Yhp1p'),:}==1))),2)));
writetable(YHP1pRxns,'output/YHP1_targetRxn.csv','WriteRowNames',1);




