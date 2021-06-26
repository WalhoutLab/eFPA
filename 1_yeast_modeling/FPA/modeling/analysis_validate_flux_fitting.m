addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/

%% 1. load the model and prepare the model
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
%% load flux
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = fluxTbl.([conditions{i},'_QP']);
end
rxnLabel = fluxTbl.Model_Reaction_ID;
dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
fluxMat_raw = fluxMat;
factor = repmat(dwTbl.gDCW_mL',size(fluxMat_raw,1),1);
fluxMat_raw = fluxMat_raw * 1000 ./ factor; %mmoles/hr/gDW
%% MOMA 
% v_ref = zeros(length(model.rxns),1);
% [A B] = ismember(model.rxns,rxnLabel);
% v_ref(A) = fluxMat_raw(B(A),1);
[solutionDel, solStatus] = linearMOMA_special_xl(model, fluxMat_raw(:,10),rxnLabel);
solutionDel.f

%% check diff
flux_QP = fluxMat_raw(:,10);
[A B] = ismember(rxnLabel,model.rxns);
flux_fit = solutionDel.full(B(A));
delFlux = abs(flux_QP - flux_fit);
checkTbl = table(rxnLabel, printRxnFormula_XL(model,rxnLabel,0),flux_QP, flux_fit, delFlux);
checkTbl = sortrows(checkTbl,5,'descend');

%% check delta flux for all 
for i = 1:size(fluxMat_raw,2)
    [solutionDel, solStatus] = linearMOMA_special_xl(model, fluxMat_raw(:,i),rxnLabel);
    deltaFlux(i) = solutionDel.f;
end
histogram(deltaFlux);
histogram(deltaFlux ./ sum(abs(fluxMat_raw),1) .* 100); % percetage change in flux
xlabel('flux diff (%)')
%% predict glucose exchange 
checkset = {'r_1166','r_1115','r_1244'};
for i = 1:size(fluxMat_raw,2)
    noGluInd = ~ismember(rxnLabel,checkset);
    [solutionDel, solStatus] = linearMOMA_special_xl(model, fluxMat_raw(noGluInd,i),rxnLabel(noGluInd));
    [A B] = ismember(checkset, model.rxns);
    pred = solutionDel.full(B(A));
    [A B] = ismember(checkset, rxnLabel);
    ref = fluxMat_raw(B(A),i);
    deltaPredPerc(:,i) = abs(pred - ref) ./ ref;
end


plot(1:25,deltaPredPerc(1,:) .*100);
ylabel('error (%)');
title('glucose transport');

plot(1:25,deltaPredPerc(2,:) .*100);
ylabel('error (%)');
title('ammonium transport');

plot(1:25,deltaPredPerc(3,:) .*100);
ylabel('error (%)');
title('phosphate transport');





