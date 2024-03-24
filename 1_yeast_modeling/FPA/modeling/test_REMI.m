function [n_corr_REMI] = test_REMI(cutoffs)
%%
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/
addpath('./../scripts/')
addpath('./../../../PlotPub/lib/')
addpath(genpath('./../other_methods/REMI/'));
% initCobraToolbox
%% model setup

% generally consistent with eFPA setup

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

% load data etc
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

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

%% make the REMI matrix 

% We run REMI with its own GRP parsing algorithm and just supply the
% fold-changes
ExpTable = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)

% we optimize the REMI parameter for fold-change thresholding in case it is
% important 
% cutoff    ncorr
% 1         17
% 1.1       28
% 1.2       47
% 1.3       32
% 1.5       25
% 1.6       2
% 1.7       3
% 1.8       8
% 1.9       6
% 2         2

% cutoffs = [1.2 1.3 1.5 1.6];
n_corr_REMI = [];
model_ori = model;
for zz = 1:length(cutoffs)
    cutoff = cutoffs(zz)
    model = model_ori;

    up_cutoff=cutoff; % more than 2 fold up regulated is known as up regulated genes (This paramter can change)
    down_cutoff=1/cutoff;% less than 1/2 fold  is known as down regulated genes 
    
    % now we are build gex Model:  which integeartes relative gene expression
    % we directly use the FC against meta-reference - so we get a matrix of
    % FC in flux against meta-reference (no problem for correlation)
    REMImat = zeros(length(model.rxns),length(conditions));

    % impose biomass as standard REMI
    model = changeRxnBounds(model,'r_2111',0.1, 'l');
    
    for index = 1:length(conditions)
        name = conditions{index};
        myExp = ExpTable.(name);
        fc = 2.^myExp;
        [rxns,ratio]=evaluateGPR(model,ExpTable.Gene,fc,@geomean,@mean);
        
        % find up- and down- regulated reactions
        indUP=find(ratio>up_cutoff);
        ratioUP=ratio(indUP);
        indDOWN=find(ratio<down_cutoff);
        ratioDOWN=ratio(indDOWN);
        regRxnRatio=[ratioUP;ratioDOWN];
        regRxns=[rxns(indUP)';rxns(indDOWN)'];
    
        % avoid numerical error (more than 100 fold is taken only 100)
        regRxnRatio(regRxnRatio>100)=100;
        regRxnRatio(regRxnRatio<1/100)=1/100;
    
        % if we want to add relative constraint into TFA model then we need to add
        % net flux variable to the model using 'addNetFluxVariablesNEW'
        % and if one want to add into FBA model  then evalute scripts given below
        mTmp=model;
        mTmp.rev=ones(numel(mTmp.rxns),1);
        use_m=addUseVariablesDH(mTmp);
        netM=addNetFluxVariablesNEW(use_m);
        [gex]=addRelConsExpression(netM,netM,regRxns,regRxnRatio);
        sol_gex=solveTFBAmodel(gex,false,'gurobi',true);
    
        remiOut = sol_gex;
        [A B] = ismember(cellfun(@(x) ['PERTURB_NF_',x],model.rxns,'UniformOutput',false),gex.varNames);
        [C D] = ismember(cellfun(@(x) ['NF_',x],model.rxns,'UniformOutput',false),gex.varNames);
        REMImat(1:length(model.rxns),index) = abs(remiOut.x.x(B(A))) ./ abs(remiOut.x.x(D(C)));
%         [A B] = ismember(cellfun(@(x) ['PERTURB_R_',x],model.rxns,'UniformOutput',false),gex.varNames);
%         [C D] = ismember(cellfun(@(x) ['R_',x],model.rxns,'UniformOutput',false),gex.varNames);
%         REMImat((1+length(model.rxns)):2*length(model.rxns),index) = remiOut.x.x(B(A)) ./ abs(remiOut.x.x(D(C)));
    end
    REMImat(isnan(REMImat)) = 1; % nan for no change
    REMImat(isinf(REMImat)) = 100; % max change is 100 fold
    REMImat(REMImat==0) = 1e-2; % min change is 0.01 fold
    
    %% correlation
    fluxMat = fluxMat_normalized;
    %Computing the correlation between a reaction expression and measured growth rate
    rho=[];
    p=[];
    testedRxn = {};
    rxnLabel = fluxTbl.Model_Reaction_ID;
    valid_rxns = model.rxns;
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            prediction = REMImat(strcmp(valid_rxns,rxnLabel{j}),:)';
            [rho(end+1),p(end+1)] = corr(prediction,abs(fluxMeasure)','type','Pearson');
        end
    end
    testedRxn_pro = testedRxn;
    p_pro = p;
    rho_pro = rho;
    figure
    histogram(rho)
    sum(rho(p<0.05)>0);
    
    fdr = mafdr(p,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    n_corr_REMI(zz) = sum(fdr<0.05 & rho > 0);
end
