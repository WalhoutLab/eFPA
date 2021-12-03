% we are preparing a final table where for each reaction, we annotate its
% basic information, its class (expression measured? expression change? expression
% qualitative? flux change? flux qualitative? predicted? predicted with
% some distance boundary?) and p, r and rho, etc. and we will make pie-chart sumary of this; with
% that, we should label/distinguish the non-correlated reaction from
% unchanged reactions.

%% first make the master information table
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

extRxns = model.rxns(findExcRxns(model));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

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
%% predictions and correlations
normalizedLevel = normalizedLevel_pro_perPro;
valid_rxns = valid_rxns_pro_perPro ;
%Computing the correlation between a reaction expression and measured growth rate
r_exp_only=[];
p_exp_only = [];
testedRxn_exp_only = {};
rxnLabel = fluxTbl.Model_Reaction_ID;
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn_exp_only(end+1) = rxnLabel(j);
        % [rho(end+1),p_rho(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
        [r_exp_only(end+1),p_exp_only(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        % 
        %Correcting for multiple hypothesis using FDR and significance
        %level of?
    end
end
fdr_r_exp = mafdr(p_exp_only,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r_exp_only(fdr_r_exp<0.05)>0));
sig_rxn_exp = testedRxn_exp_only(r_exp_only>0 & fdr_r_exp<0.05);

% correlation of FPA prediction
load(['output/Titration_relativeExp_wtdDist_expDecay.mat'])
nn = 9;
zz = 2;
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
r_FPA=[];
deltaminmax = [];
testedRxn_FPA = {};
p_FPA = [];
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn_FPA(end+1) = rxnLabel(j);
        [r_FPA(end+1),p_FPA(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end
            
fdr_r_FPA = mafdr(p_FPA,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r_FPA(fdr_r_FPA<0.05)>0));
sig_rxn_FPA = testedRxn_FPA(r_FPA>0 & fdr_r_FPA<0.05 & deltaminmax > 0.2);


%% supp1 A: reaction information
% reaction ID, formula, GPR, expression change/measured?, flux
% change/measured?, correlated?, predicted? optimal FPA r? optimal FPA FDR?
rxnTbl = struct();
rxnTbl.rxnID = rxnLabel;
rxnTbl.formula = printRxnFormula_XL(model,rxnLabel,0);
rxnTbl.GPR = printGPRForRxns(model,rxnLabel,0);
rxnTbl.expression_type = repmat({'Not Assigned'},length(rxnLabel),1);
rxnTbl.flux_type = repmat({'Not Assigned'},length(rxnLabel),1);
rxnTbl.correlated = repmat({'Not Assigned'},length(rxnLabel),1);
rxnTbl.predicted_by_default_FPA = repmat({'Not Assigned'},length(rxnLabel),1);
rxnTbl.PCC_by_target_expression = nan(length(rxnLabel),1);
rxnTbl.FDR_by_target_expression = nan(length(rxnLabel),1);
rxnTbl.PCC_by_default_FPA = nan(length(rxnLabel),1);
rxnTbl.FDR_by_default_FPA = nan(length(rxnLabel),1);
rxnTbl.predicted_by_FPA = repmat({'Not Assigned'},length(rxnLabel),1);
heatmapTbl = readtable('output/relCorr_heatmapTbl_realDist.csv');


for i = 1:length(rxnLabel)
    myrxn = rxnLabel{i};
    if ~any(strcmp(myrxn,valid_rxns_pro_perPro))
        if isempty(rxnTbl.GPR{i}) 
            rxnTbl.expression_type{i} = 'No GPR';
        else
            rxnTbl.expression_type{i} = 'Not measured';
        end
    else
        myexp = normalizedLevel_pro_perPro(strcmp(myrxn,valid_rxns_pro_perPro),:);
        rxnTbl.expression_type{i} = checkVariability(myexp);
    end
   
    myflux = fluxMat_normalized(i,:);
    rxnTbl.flux_type{i} = checkVariability(myflux);
    if any(strcmp(sig_rxn_exp,myrxn))
        rxnTbl.correlated{i} = 'Yes';
    else
        rxnTbl.correlated{i} = 'No';
    end
    if any(strcmp(sig_rxn_FPA,myrxn))
        rxnTbl.predicted_by_default_FPA{i} = 'Yes';
    else
        rxnTbl.predicted_by_default_FPA{i} = 'No';
    end
    if any(strcmp(testedRxn_exp_only,myrxn))
        rxnTbl.PCC_by_target_expression(i) = r_exp_only(strcmp(testedRxn_exp_only,myrxn));
        rxnTbl.FDR_by_target_expression(i) = fdr_r_exp(strcmp(testedRxn_exp_only,myrxn));
    else
        rxnTbl.PCC_by_target_expression(i) = NaN;
        rxnTbl.FDR_by_target_expression(i) = NaN;
    end
    if any(strcmp(testedRxn_FPA,myrxn))
        rxnTbl.PCC_by_default_FPA(i) =  r_FPA(strcmp(testedRxn_FPA,myrxn));
        rxnTbl.FDR_by_default_FPA(i) = fdr_r_FPA(strcmp(testedRxn_FPA,myrxn));
    else
        rxnTbl.PCC_by_default_FPA(i) = NaN;
        rxnTbl.FDR_by_default_FPA(i) = NaN;
    end
    
    if any(strcmp(heatmapTbl.Row,myrxn))
        rxnTbl.predicted_by_FPA{i} = 'Yes';
    else
        rxnTbl.predicted_by_FPA{i} = 'No';
    end
     
end

rxnTbl = struct2table(rxnTbl);
writetable(rxnTbl,'output/summary_table_reaction_information.csv');

% enrichment for high diversity flux 
tabulate(rxnTbl.expression_type(strcmp(rxnTbl.predicted_by_FPA,'No')))
tabulate(rxnTbl.expression_type(strcmp(rxnTbl.predicted_by_FPA,'Yes')))

%% supp1 B: flux information
mole_per_mlcell_perHour = array2table(fluxMat);
mole_per_mlcell_perHour.Properties.VariableNames = conditions;
mole_per_mlcell_perHour.Properties.RowNames = rxnLabel;
writetable(mole_per_mlcell_perHour,'output/supp1B_mole_per_mlcell_perHour.csv','WriteRowNames',true);

normalizedFlux = array2table(fluxMat_normalized);
normalizedFlux.Properties.VariableNames = conditions;
normalizedFlux.Properties.RowNames = rxnLabel;
writetable(normalizedFlux,'output/supp1B_normalizedFlux.csv','WriteRowNames',true);

rawFlux = array2table(fluxMat_raw);
rawFlux.Properties.VariableNames = conditions;
rawFlux.Properties.RowNames = rxnLabel;
writetable(rawFlux,'output/supp1B_rawFlux.csv','WriteRowNames',true);
%% supp1 C: original expression (copy from the SIMMER TABLE './../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx')
%% supp1 D: normalized reaction expression
normalizedExp = array2table(normalizedLevel_pro_perPro);
normalizedExp.Properties.VariableNames = conditions;
normalizedExp.Properties.RowNames = valid_rxns_pro_perPro;
writetable(normalizedExp,'output/supp1D_normalizedExpression.csv','WriteRowNames',true);

%% supp1 E: relative flux potentials (base2 d6)
relFluxPotential = array2table(relFP);
relFluxPotential.Properties.VariableNames = conditions;
relFluxPotential.Properties.RowNames = valid_rxns;
writetable(relFluxPotential,'output/supp1E_relFluxPotential.csv','WriteRowNames',true);

%% supp1 F: PCCs - skip for now
%% supp1 G: p-values for the PCCs - skip for now
%%
% this table will also tell how many reactions are predicted to change but
% not change/correct; so this is a estimation of false positive by the
% local pathway rule ==> we dont analyze this for now

% how many reactions are changing significantly in expression but not
% coherent with flux
%% this is an important result! High diversity expression has a better correlation with flux 
tabulate(rxnTbl.correlated)
tabulate(rxnTbl.correlated(strcmp(rxnTbl.expression_type,'High Diversity')))
tabulate(rxnTbl.predicted_by_default_FPA(strcmp(rxnTbl.expression_type,'High Diversity')))
tabulate(rxnTbl.predicted_by_FPA(strcmp(rxnTbl.expression_type,'High Diversity')))
%% how about high diversity flux? it is also true! so the current number is still a underestimation of how much expression can proxy flux
tabulate(rxnTbl.correlated)
tabulate(rxnTbl.correlated(strcmp(rxnTbl.flux_type,'High Diversity')))
tabulate(rxnTbl.predicted_by_default_FPA(strcmp(rxnTbl.flux_type,'High Diversity')))
tabulate(rxnTbl.predicted_by_FPA(strcmp(rxnTbl.flux_type,'High Diversity')))
%%
% for reactions with expression and flux changes but not predicted any,
% check for if pathway level is consistent but not correlating with flux; 

% we also look into features of regulated and not regulated rxns ((1)not changes 
% (2) not correlated); for
% example, is all not regulated rxns are genes that are associated with
% multiple pathways?