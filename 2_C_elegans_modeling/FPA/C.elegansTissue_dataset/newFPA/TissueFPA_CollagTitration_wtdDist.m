%% This is to reproduce the collagen FPA titration in Fig.5D
% THIS SCRIPT IS NOT WELL-COMMENTED AND IS ONLY FOR REPRODUCING RESULTS

%% env variables
addpath /share/pkg/gurobi/900/linux64/matlab/
addpath ~/cobratoolbox/
addpath ./input/
addpath ./../../scripts/
addpath ./../../scripts/oriMERGE/
addpath ./../../../input/GEMs/
addpath ./../../../bins/
initCobraToolbox(false)
%% prepare files
load('Tissue.mat');
% the default constriants for dual model
model = changeRxnBounds(model,'EXC0050_L',-1000,'l');
model = changeRxnBounds(model,'EX00001_E',-1000,'l');
model = changeRxnBounds(model,'EX00007_E',-1000,'l');
model = changeRxnBounds(model,'EX00009_E',-1000,'l');
model = changeRxnBounds(model,'EX00080_E',-1000,'l');
model = changeRxnBounds(model,'RM00112_I',0,'l');
model = changeRxnBounds(model,'RM00112_X',0,'l');
model = changeRxnBounds(model,'RM00112_I',0,'u');
model = changeRxnBounds(model,'RM00112_X',0,'u');
model = changeRxnBounds(model,'DMN0033_I',0,'l');
model = changeRxnBounds(model,'DMN0033_X',0,'l');
model = changeRxnBounds(model,'DMN0033_I',1,'u');
model = changeRxnBounds(model,'DMN0033_X',1,'u');
model = changeRxnBounds(model,'RMC0005_I',0,'l');
model = changeRxnBounds(model,'RMC0005_X',0,'l');
model = changeRxnBounds(model,'RMC0005_I',0,'u');
model = changeRxnBounds(model,'RMC0005_X',0,'u');
% parseGPR takes hugh amount of time, so preparse and integrate with model
% here 
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;

% this loads the raw directional distance matrix, first need to take the
% minimum for future use
% to load from output of distance calculator, uncomment the following line
% distance_raw = readtable('distance_fromRxnToRxn.tsv','FileType','text','ReadRowNames',true);
distance_raw = readtable('input/distanceMatrix_weighted.txt','FileType','text','ReadRowNames',true);
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

% set the inter-tissue distance to zero
X_ind = contains(labels,'_X_');
I_ind = contains(labels,'_I_');
distMat(X_ind,I_ind) = 0;
distMat(I_ind,X_ind) = 0;

% load the special penalties
manualPenalty = table2cell(readtable('manualPenalty.csv','ReadVariableNames',false,'Delimiter',','));
% load the special distance
manualDist = table2cell(readtable('manualDist.csv','ReadVariableNames',false,'Delimiter',','));

% load the confidence level
confidenceTbl = readtable('confidenceTable.tsv','FileType','text','ReadRowNames',true);
confidenceTbl.Properties.RowNames = cellfun(@(x) [x(1:end-1),'_',x(end)],confidenceTbl.Properties.RowNames,'UniformOutput',false);

%% prepare expression files
expressionTbl = readtable('expressionTable.tsv','FileType','text','ReadRowNames',true);
tissueLabel = expressionTbl.Properties.VariableNames;
master_expression = {};
for i = 1:length(tissueLabel)
    expression = struct();
    % NOTE: the expression (penalty) of _I compartment rxns is handled
    % specially by the constantPenalty
    expression.genes = cellfun(@(x) [x,'_X'], expressionTbl.Properties.RowNames, 'UniformOutput', false);
    expression.value = expressionTbl.(tissueLabel{i});
    master_expression{i} = expression;
end
%% constrain model
model = changeRxnBounds(model,'RM00112_X',1000,'u');
model = changeRxnBounds(model,'RM00112_I',1000,'u');
model = changeRxnBounds(model,'RM00112_X',-1000,'l');
model = changeRxnBounds(model,'RM00112_I',-1000,'l');

bacMW=966.28583751;
model.S(end, strcmp('EXC0050_L',model.rxns)) = 0.02*bacMW*0.01; %the side cap
model.S(ismember(model.mets,{'atp_I[c]','h2o_I[c]','adp_I[c]','h_I[c]','pi_I[c]'}), strcmp('DGR0007_L',model.rxns)) = 0; 

%remove the NGAM 
model = changeRxnBounds(model,'RCC0005_X',0,'l');
model = changeRxnBounds(model,'RCC0005_I',0,'l');
%% block the -2/-3 reactions
model_irrev = convertToIrreversible(model);
%unify naminclature
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); %has no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
X_rxns = model_irrev.rxns(contains(model_irrev.rxns,'_X_'));
for i = 1:length(tissueLabel)
    blockList{i} = confidenceTbl.Properties.RowNames((confidenceTbl.(tissueLabel{i}) <= -2));
    if strcmp(tissueLabel{i},'Intestine')
        blockList{i} = [blockList{i};X_rxns]; %block all X tissue reactions for intestine
    end     
end
blockList{end+1} = {};%for supercond, dont block anything when optimizing X tissue
blockList_I = blockList;
blockList_I{end} = X_rxns;%still block all X tissue when optimzing supercond of a Intestine tissue
%% prepare paremeters
nSeq = 0:0.5:10; % the distance order sequence to titrate
maxDist = max(distMat(~isinf(distMat)));

changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

myCluster = parcluster('local');
myCluster.NumWorkers = 128;
saveProfile(myCluster);
%parpool(4,'SpmdEnabled',false);
%% setup special penalties for X tissue and I tissue
manualPenalty2 = table2cell(readtable('nonIntestineManualPenalty.csv','ReadVariableNames',false,'Delimiter',','));
manualPenalty2 = [manualPenalty2;manualPenalty];
master_expression2 = {};
for i = 1:length(tissueLabel)
    expression = struct();
    % NOTE: the expression (penalty) of _I compartment rxns is handled specially by the constantPenalty
    expression.genes = cellfun(@(x) [x,'_I'], expressionTbl.Properties.RowNames, 'UniformOutput', false);
    expression.value = expressionTbl.(tissueLabel{i});
    master_expression2{i} = expression;
end
penalty = calculatePenalty_oriMERGE(model,master_expression2,manualPenalty2);
I_ind = contains(model.rxns,'_I');
constantPenalty2 = [model.rxns(I_ind),mat2cell(penalty(I_ind,strcmp(tissueLabel,'Intestine')),ones(sum(I_ind),1))];
penalty = calculatePenalty_oriMERGE(model,master_expression2,manualPenalty);
I_ind = contains(model.rxns,'_I');
constantPenalty = [model.rxns(I_ind),mat2cell(penalty(I_ind,strcmp(tissueLabel,'Intestine')),ones(sum(I_ind),1))];

%% titrate the Flux Potentials for collagen
targetExRxns = {'NewMet_collg[c]'}; %it is metabolite-centered FPA
% the difference is specific rxn(s) need to be blocked to prevent loop and
% false FPA value. these specific rxns are related to the objective in the
% optimization (producing or degradating); hence, we need to optimize the
% forward and/or reverse direction of a transporter seperately.
% we are mainly interested in the transportation of organic molecules. some
% co-transporting metabolite are not the target of interest so that they
% shouldnt be blocked(limiting)
helperMet = {'co2','amp','nadp','nadph','ppi','o2','nadh','nad','pi','adp','coa', 'atp', 'h2o', 'h', 'gtp', 'gdp','etfrd','etfox','na1','oh1','cl','i'}; %define a set of co-transported metabolite that are not the metabolite of interest
penalty_defined_X = calculatePenalty_oriMERGE(model,master_expression,manualPenalty2);
penalty_defined_I = calculatePenalty_oriMERGE(model,master_expression,manualPenalty);
fluxEfficiency_X_met = cell(length(nSeq),size(penalty,2));
fluxEfficiency_plus_X_met = cell(length(nSeq),size(penalty,2));
fluxEfficiency_I_met = cell(length(nSeq),size(penalty,2));
fluxEfficiency_plus_I_met = cell(length(nSeq),size(penalty,2));
environment = getEnvironment();
parfor i = 1:length(nSeq)
    restoreEnvironment(environment);
    n = nSeq(i);
     % it is like a demand reaction except for adding a new reaction
    myMet = regexprep(targetExRxns{1},'^NewMet_','');
    % find the associated end reactions (producing the metabolite)
    myRxns = model.rxns(any(model.S(cellfun(@(x) ~isempty(regexp(x,['^',myMet(1:end-3),'_.\[.\]$'],'once')),model.mets),:),1)); %all rxns use these mets
    myRxns = myRxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE)','once')),myRxns));
    % block these reactions and optimize 
    % only block the reactions for corresponding tissue
    model_tmp = model;
    % add the new demand 
    model_tmp = addReaction(model_tmp,['DMN',myMet,'_X'],'reactionFormula',[myMet(1:end-3),'_X',myMet(end-2:end),' ->'],'geneRule', 'NA','printLevel',0);
    targetrxn_fullName = ['DMN',myMet,'_X'];
    % the distance of this new demand rxn needs to be calculated; Note
    % that the distance is actually 1+min(distance(any rxn that
    % contains this met)
    [distMat_raw_2,labels_X_2] = updateDistance(targetrxn_fullName,model_tmp,distMat_raw,labels,maxDist);
    % take the min and block intestine
    distMat_X_2 = distMat_raw_2;
    for ii = 1:size(distMat_X_2,1)
        for jj = 1:size(distMat_X_2,2)
            distMat_X_2(ii,jj) = min([distMat_raw_2(ii,jj),distMat_raw_2(jj,ii)]);
        end
    end
    % set the inter-tissue distance to zero
    X_ind = contains(labels_X_2,'_X_');
    I_ind = contains(labels_X_2,'_I_');
    distMat_X_2(X_ind,I_ind) = 0;
    distMat_X_2(I_ind,X_ind) = 0;

    penalty_defined_X_2 = [penalty_defined_X; ones(1,size(penalty_defined_X,2))];
    myRxns_X = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_X)$','once')),myRxns));
    model_tmp.ub(ismember(model_tmp.rxns,myRxns_X)) = 0;
    model_tmp.lb(ismember(model_tmp.rxns,myRxns_X)) = 0;
    [fluxEfficiency_X_met(i,:),fluxEfficiency_plus_X_met(i,:)] = FPA_oriMERGE(model_tmp,{targetrxn_fullName},master_expression,distMat_X_2,labels_X_2,n, manualPenalty2,manualDist,maxDist,blockList,constantPenalty2,false,penalty_defined_X_2);

    model_tmp = model;
    model_tmp = addReaction(model_tmp,['DMN',myMet,'_I'],'reactionFormula',[myMet(1:end-3),'_I',myMet(end-2:end),' ->'],'geneRule', 'NA','printLevel',0);
    targetrxn_fullName = ['DMN',myMet,'_I'];
    [distMat_raw_2,labels_I_2] = updateDistance(targetrxn_fullName,model_tmp,distMat_raw,labels,maxDist);
    % take the min and block intestine
    distMat_I_2 = distMat_raw_2;
    for ii = 1:size(distMat_I_2,1)
        for jj = 1:size(distMat_I_2,2)
            distMat_I_2(ii,jj) = min([distMat_raw_2(ii,jj),distMat_raw_2(jj,ii)]);
        end
    end
    % set the inter-tissue distance to zero
    X_ind = contains(labels_I_2,'_X_');
    I_ind = contains(labels_I_2,'_I_');
    distMat_I_2(X_ind,I_ind) = 0;
    distMat_I_2(I_ind,X_ind) = 0;
    penalty_defined_I_2 = [penalty_defined_I; ones(1,size(penalty_defined_I,2))];
    myRxns_I = myRxns(cellfun(@(x) ~isempty(regexp(x,'(_L|_I)$','once')),myRxns));
    model_tmp.ub(ismember(model_tmp.rxns,myRxns_I)) = 0;
    model_tmp.lb(ismember(model_tmp.rxns,myRxns_I)) = 0;
    [fluxEfficiency_I_met(i,:),fluxEfficiency_plus_I_met(i,:)] = FPA_oriMERGE(model_tmp,{targetrxn_fullName},master_expression,distMat_I_2,labels_I_2,n, manualPenalty,manualDist,maxDist,blockList_I,constantPenalty,false,penalty_defined_I_2);
end
%% get the relative FP
fluxEfficiency_X = [fluxEfficiency_X_met];
fluxEfficiency_I = [fluxEfficiency_I_met];
% get the final result of reltaive flux potential
relFP_f = nan(size(fluxEfficiency_X,1),length(tissueLabel));%flux potential for forward rxns
relFP_r = nan(size(fluxEfficiency_X,1),length(tissueLabel));%flux potential for reverse rxns
for i = 1:size(fluxEfficiency_X,1)
    for j = 1:length(tissueLabel)
        if ~strcmp(tissueLabel{j},'Intestine')
            relFP_f(i,j) = fluxEfficiency_X{i,j}(1) ./ fluxEfficiency_X{i,end}(1);
            relFP_r(i,j) = fluxEfficiency_X{i,j}(2) ./ fluxEfficiency_X{i,end}(2);
        else
            relFP_f(i,j) = fluxEfficiency_I{i,j}(1) ./ fluxEfficiency_I{i,end}(1);
            relFP_r(i,j) = fluxEfficiency_I{i,j}(2) ./ fluxEfficiency_I{i,end}(2);
        end
    end
end
% save('FPA_rxns_collag.mat','relFP_f','relFP_r');
% save('FPA_workspace_3.mat');
%% plot and evaluate (Fig 5D)
figure(1)
hold on
for i = 1:length(tissueLabel)
    plot(nSeq, relFP_f(:,i),'.-','LineWidth',5);
end
legend(tissueLabel);
    
    