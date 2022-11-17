%% README FOR REVIEWERS
% Dear reviewers, this walkthrough is a temporary replacement for the final 
% user tutorial that we are developing. We aim to provide a more user-
% friendly eFPA/FPA master function with application demos to several common 
% network models (e.g. allowing for flexible input/output, easy tuning of 
% parameters, one-click solution of metabolite-centric eFPA, etc). The 
% user-friendly eFPA package will be completed in a few months. This 
% temporary demo is only for a test of the algorithm and only contains the 
% demonstration of the most basic use of eFPA. 

%% Add required path 
% **** make sure you do this every time you run the demo!****
% add path for required functions/inputs
addpath ./input/
addpath ./functions/
addpath ./input/
initCobraToolbox(false);


%% DEMO: APPLICATION TO THE GENERIC C. ELEGANS MODEL WITH BULK RNA-SEQ INPUT
%% 1. Prepare the model
load('iCEL1314.mat');
% Users may add their own constraints here (i.e., nutritional input constraints)
% the following constraints are specific for C. elegans simulations
model = changeRxnBounds(model,'EXC0050',-1000,'l');% we allow unlimited bacteria uptake
model = changeRxnBounds(model,'RCC0005',0,'l');% remove the NGAM 
model.S(ismember(model.mets,{'atp[c]','h2o[c]','adp[c]','h[c]','pi[c]'}), strcmp('DGR0007',model.rxns)) = 0;% remove the energy cost for bacteria digestion

%% 2. Load the expression files, distance matrix, and other optional inputs

% ```load expression files```
% Expression matrix can be in plain text and in any normalized
% quantification metric like TPM or FPKM.
% We use the RNA-seq data from Bulcha et al, Cell Rep (2019) as an example
expTbl = readtable('exampleExpression.csv');
% For demo purpose, we only analyze the FPA of four conditions in the
% expression dataset.
conditions = {'N2_OP50', 'N2_B12', 'nhr10_OP50','nhr10_B12'};

% ```preprocess the expression table```
% Since expression tables from different data source may have different
% format, we require users to re-organize their table according to the
% following procedures.
% In eFPA, we convert the expression table into a variable called
% "master_expression".
% ```REQUIREMENT OF "master_expression"```
% To standardize the format and symbols of different expression matrix, we 
% use a "master_expression" variable as the expression input of eFPA. The
% format requirements are as follows:
% (1) it must be a cell arrary of structure variables;
% (2) each structure variable contains the expression information of one
%     condition; and the order of the structure variables in the cell array
%     determines the order of FP values in the output;
% (3) each structure variable must contain two fields: "genes" and "value".
%     "genes" and "value" should be of equal length.
% (4) every input "genes" should be measured in all conditions. In other
%     words, it is not recommend to have different "genes" list in
%     different conditions (although it is allowed by eFPA).
%
% Make a new master_expression for these four conditions.
master_expression = {};
% Get the index of genes in the model
geneInd = ismember(expTbl.Gene_name, model.genes); 
for i = 1:length(conditions)
    expression = struct();
    expression.genes = expTbl.Gene_name(geneInd);
    expression.value = expTbl.(conditions{i})(geneInd);
    master_expression{i} = expression;
end

% ```load the distance matrix```
% We load from the output of the distance calculator. 
% For usage of distance calculator, please refer to the `3_distance_calculation` folder
distance_raw = readtable('distanceMatrix_weighted.txt','FileType','text','ReadRowNames',true); 
% convert the directional distance to undirection distance by taking the
% minimal distance
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

% ```load the special penalties```
% - We set penalty for all Exchange, Demand, Degradation and Sink reactions
%  to 0, to not penalize the external reactions
manualPenalty = table2cell(readtable('manualPenalty_generic.csv','ReadVariableNames',false,'Delimiter',','));

%% 3. Run regular FPA analysis

% The FPA/eFPA is designed with parfor loops for better speed, so we first 
% initialize the parpool
parpool(4)

% ```setup some basic parameters for eFPA```
distBound = 6; % distance boundary; in this demo, we only used the exponential decay formula; for binary decay, refer to the yeast section or wait for the final tutorial
changeCobraSolverParams('LP','optTol', 10e-9); % solver parameter
changeCobraSolverParams('LP','feasTol', 10e-9); % solver parameter

% ```setup target reactions```
% we perform eFPA analysis for three reactions (PP shunt pathway) as an example
targetRxns = {'RM04432';'RM03045'};

% To run eFPA, simply write:
[FP,~] = eFPA(model,targetRxns,master_expression,distMat,labels,distBound, manualPenalty);
%
% Calculate relative flux potential (rFP) 
relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end

%% 4. inspect the eFPA result
figure(1)
c = categorical(regexprep(conditions,'_','-'));
bar(c,relFP_f(1,:))
title('rFP of Propanoyl-CoA:(acceptor) 2,3-oxidoreductase flux')

figure(2)
c = categorical(regexprep(conditions,'_','-'));
bar(c,relFP_f(2,:))
title('rFP of 3-Hydroxypropionyl-CoA hydrolyase flux')

% We can see that the eFPA prediction recaptures the repressing of
% Propanoyl-CoA:(acceptor) 2,3-oxidoreductase flux by vitamin B12
% treatment, and the loss of flux activation after nhr10 is deleted.
% (Bulcha et al, Cell Rep (2019), PMID: 30625328).

