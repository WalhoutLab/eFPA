%% This is to reproduce the C. elegans tissue FPA results in supp table S7 and S8
%% set up the env variables
addpath ~/cobratoolbox/
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
initCobraToolbox(false)
%% prepare model
load('input/ihuman_COBRA.mat');
model.subSystems = [model.subSystems{:}]';
% the default constraints are unlimited, so we dont change anything

% ```special treatments for iHuman1```
% Remove parentathsis in the reaction ID (which will be changed in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% change the name starting with numbers 
model.rxns = regexprep(model.rxns,'(^[0-9])','x$1');


% find the target internal reactions 
exludeSys = {'Transport reactions','Exchange/demand reactions'};
intRxns = model.rxns(~ismember(model.subSystems,exludeSys));

% create missing field
model = creategrRulesField(model);
%%
parpool(4);
%%
environment = getEnvironment();
parfor i = 1:length(model.rxns)
    restoreEnvironment(environment);
    model1 = changeObjective(model,model.rxns(i));
    tmp = optimizeCbModel(model1,'max');
    UBs(i) = tmp.obj;
    tmp = optimizeCbModel(model1,'min');
    LBs(i) = tmp.obj;
end
%%
lowFluxRxns1 = model.rxns(UBs > 0 & UBs <= 1);
lowFluxRxns2 = model.rxns(LBs < 0 & LBs >= -1);
lowFluxRxns = union(lowFluxRxns1,lowFluxRxns2);
save('lowFluxRxns.mat','lowFluxRxns');
