function [n_corr_deltaFBA] = test_deltaFBA()
% Set paths, initiate cobra and solver parameters
addpath(genpath('./../other_methods/DeltaFBA/'));
[cobraParams,solverParams] = parseSolverParameters('MILP');
cobraParams.timeLimit = 1000;
cobraParams.printLevel = 1;

% Load model and prepare constraints here 
% constraints should be well set up to perform regular FBA
model = loadYeatModel();
% % the following nutrients need to be set manually
% % to start with basic FPA, we allow unlimited exchange
% % phosphate exchange
model.lb(strcmp(model.rxns,'r_2005')) = -1000;
% glucose exchange
model.lb(strcmp(model.rxns,'r_1714')) = -1000;
% ammonium exchange 
model.lb(strcmp(model.rxns,'r_1654')) = -1000;
% uracil 
model.lb(strcmp(model.rxns,'r_2090')) = -1000;
% leucine
model.lb(strcmp(model.rxns,'r_1899')) = -1000;
% maintanence - remove maintanance to be consistent with the similar case in the original paper (Tienen GSE19420)
model = changeRxnBounds(model,'r_4046',0,'l'); % maintance 
model = changeRxnBounds(model,'r_4046',1000,'u'); % maintance 

% Read model and optimze Wild type model
model = generateRules(model);
FBA = optimizeCbModel(model);

% no imposed biomass as we skip the exp delta flux constraints
% it is unclear how to define a media constraint set to fit various
% limiting conditions, so we simply go with imposing biomass production
% model = changeRxnBounds(model,'r_2111',0.1, 'l');

% Find minimized model flux results to Set maximum flux of a single reaction

[MinimizedFlux modelIrrev]= minimizeModelFlux(model);

% add the irreversible flag
revRxns = setdiff(modelIrrev.rxns, model.rxns);
uniqueRxns = setdiff(unique(regexprep(revRxns,'_(b|f)$','')), 'netFlux');
% setdiff(strcat(uniqueRxns,'_f'), revRxns)
model.rev = ismember(model.rxns, uniqueRxns);

minflux = optimizeCbModel(modelIrrev,'min'); % should be zero as there is no maintanance and biomass constriant

if isfield(minflux, 'full')==1
    vec = minflux.full;
else isfield(minflux,'x')== 1
    vec = minflux.x;
end

% check for flux scale under a free exchange constraint
maxflux_val = 200; % max(abs(vec(1:(end-1)))); % this is a free constrain model, so we use the parameter used in the similar case in the original paper (Tienen GSE19420)

% Check for essential reactions that are constrained
essential_idx = union(find(model.ub~=0 & model.ub<1000), find(model.lb~=0 & model.lb>-1000));
% Note = Check essential reactions and prune set

% Check for reactions that have no flux flows
zeroflow_idx = find(model.lb==0 & model.ub==0); %zero fluxes by model definition

% Create delta model by setting bounds of all reactions
[model_del, nochange_idx] = createDeltaModel(model, [], zeroflow_idx,maxflux_val);

% Set parameters for Delta FBA
epsilon = 10; % 0.1;% this is a free constrain model, so we use the parameter used in the similar case in the original paper (Tienen GSE19420)
M_prime =1e5;

% Add NetFlux Variable
[model_delNet] = addNetDeltaFlux(model_del, model.rxns);

ExpTable = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = ExpTable.Properties.VariableNames(2:end);
genes = ExpTable.Gene;
genes_data = 2.^ExpTable{:,2:end};

fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;

% flux is in  (moles / hr / mL cells); no conversion is needed. 
% in fact, correlation got worse if we normalzie the flux to / gDW first!
fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);

flux_list = rxnLabel;

P_D = zeros(numel(flux_list),numel(conditions)); %Predicted flux delta
Sol_stat = zeros(numel(conditions)-1,1); %1 is optimization for l2 norm is successful
Sol_stat2 = zeros(numel(conditions)-1,1); % MIPGAP for l2 norm

for j = 1:numel(conditions)

    %% Calculate ratio and filter
    Generatio=nan(numel(model.genes),1);
    for i=1:numel(model.genes)
        [~,ind]=ismember(model.genes{i},genes);
        if ind>0
            Generatio(i)=genes_data(ind,j);
        end
    end
    expr_genes = model.genes;
    sel = find(~isnan(Generatio));
    Generatio = Generatio(sel);
    expr_genes = expr_genes(sel);
    expressionData.gene = expr_genes;
    expressionData.value = Generatio;
    
    RxnExpr = mapExpressionToReactions(model,expressionData);
    [~,ids] = ismember(flux_list, model.rxns);
    RxnExpr_mod = RxnExpr(ids);
    RxnExpr_mod(find(RxnExpr_mod==(-1))) = 0; % from the original example code - unclear what this is for
    
    de_exprrxns = mapExpressionToReactions(modelIrrev,expressionData);
    de_rxns = modelIrrev.rxns;
    sel = find(~(de_exprrxns==(-1) | de_exprrxns==0)); % from the original example code - unclear what this is for 
    de_exprrxns = de_exprrxns(sel);
    de_rxns = de_rxns(sel);
    
    %% Define up and down regulated reactions
    cut_off = 0.25; % 0.05; % top five percentage change - as default but can be titrated 

    [~,de_indUP]=maxk(de_exprrxns,round(cut_off*size(find(de_exprrxns>1),1)));
    [~,de_indDOWN]=mink(de_exprrxns,round(cut_off*size(find(de_exprrxns<1),1)));
    
    regRxns=[de_rxns(de_indUP); de_rxns(de_indDOWN)];
    regRxnsRatio = [de_exprrxns(de_indUP) ; de_exprrxns(de_indDOWN)];
    
    %% Create binary variables for all reactions 

    % start with strigent (not quantitative) and may move to the statement
    % B to be aligned with the purpose of correlation

    % Choose 1 of the following two statements
    
    % Statement A : Relaxed conditions
    model_Bin = createBinaryUseVariable(model_delNet, epsilon, M_prime, regRxns, ones(numel(regRxns),1));
    
    % Statement B : Stringent conditions (when using this condition, please
    % make sure that generatio isnt out of bounds
    % regRxnsRatio(find(abs(regRxnsRatio)<0.1)) = 0.1;
    % regRxnsRatio(find(abs(regRxnsRatio)>100)) = 100;
    % model_Bin = createBinaryUseVariable(model_delNet, epsilon, M_prime, regRxns, regRxnsRatio);
    % 
    %% Creating objectives for maximizing the consistency of DE reactions
    % Choose 1 of the following two statements
    
    % Statement A : Without weighting the objective function with gene
    % expression changes
    [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 0);
    
    % Statement B : With weighting the objective function with gene
    % expression changes
    % de_exprrxns(find(abs(de_exprrxns)<0.1)) = 0.1;
    % de_exprrxns(find(abs(de_exprrxns)>100)) = 100;
    % [MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 1);

    [~,NFind]=ismember(strcat('NF_',model.rxns),MILPstructure.varNames);
    
    %% try add in the delta flux for exchanges and use fc agaist the first condition
    % we can then simply set the delta flux of the first condition as zero
    % in correlation

    %% Solve for 1. Maximizing consistency in DE reactions
    %Setup parameters for optimizer in gurobi
    params.OutputFlag = 1;
    params.DisplayInterval = 5;
    params.TimeLimit = cobraParams.timeLimit;
    params.MIPGap = 1e-3;
    %params.IntFeasTol = 1e-09;
    %params.FeasibilityTol = cobraParams.feasTol;
    %params.OptimalityTol = cobraParams.optTol;
    
    sol_Z_de = gurobi(MILPstructure, params);
    
    %% Solve for 2. Minimizing the 2-norm of dFBA results
    repMILP = MILPstructure;
    [~,num_vars] = size(repMILP.A);
    repMILP.lb(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
    repMILP.ub(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));

    Qmat = sparse(num_vars,num_vars);
    for i = 1:numel(NFind)
        Qmat(NFind(i), NFind(i)) = 1;
    end

    repMILP = rmfield(repMILP,'obj');
    repMILP.Q = Qmat;
    repMILP.modelsense = 'min';
    % sol_rep = gurobi(repMILP, params);

    [MILP] = create1normSolutionMILP(MILPstructure, NFind, sol_Z_de.objval);
    sol_rep = gurobi(MILP,params);

    if isfield(sol_rep,'x')
        % sol_rep = sol_rep2;
        Sol_stat(j) = 1;
        Sol_stat2(j) = sol_rep.mipgap;
    else
        sol_rep = sol_Z_de;
    end
    
     %% save results
    [~,ids] = ismember(strcat('NF_',flux_list), repMILP.varNames);
    Predicted = sol_rep.x(ids);
    P_D(:,j) = Predicted;    
end

fluxMat = fluxMat_normalized;
%Computing the correlation between a reaction expression and measured growth rate
rho=[];
p=[];
testedRxn = {};
rxnLabel = fluxTbl.Model_Reaction_ID;
valid_rxns = flux_list;
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn(end+1) = rxnLabel(j);
        prediction = P_D(strcmp(valid_rxns,rxnLabel{j}),:)';
        [rho(end+1),p(end+1)] = corr(abs(prediction),abs(fluxMeasure)','type','Pearson');
        % Correcting for multiple hypothesis using FDR and significance level of
    end
end
testedRxn_pro = testedRxn;
p_pro = p;
rho_pro = rho;
figure
histogram(rho)
sum(rho(p<0.05)>0)

fdr = mafdr(p,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
n_corr_deltaFBA = sum(fdr<0.05 & rho > 0);
