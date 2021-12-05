%% set up the env variables
addpath ~/cobratoolbox/
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
initCobraToolbox(false)
%% prepare model
load('input/ihuman_serum.mat');
model.subSystems = [model.subSystems{:}]';
% the default constraints are unlimited, so we dont change anything

% ```special treatments for iHuman1```
% Remove parentathsis in the reaction ID (which will be changed in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% change the name starting with numbers 
model.rxns = regexprep(model.rxns,'(^[0-9])','x$1');

% correct GPR
model = changeGeneAssociation(model, 'HMR_4137',...
    'ENSG00000091140 and ENSG00000110435 and (ENSG00000131828 or ENSG00000163114) and ENSG00000150768 and ENSG00000168291');


% find the target internal reactions 
exludeSys = {'Transport reactions','Exchange/demand reactions'};
intRxns = model.rxns(~ismember(model.subSystems,exludeSys));

% create missing field
model = creategrRulesField(model);
reactions = model.rxns;
environment = getEnvironment();
fprintf(['\n' repmat('.',1,length(reactions)) '\n\n']);
parfor i = 1: length(reactions)
    restoreEnvironment(environment);
    testModel = model;
    testModel = changeObjective(testModel,reactions(i));
    Vm_f = optimizeCbModel(testModel, 'max');
    FVA_f(i) = Vm_f.obj;
    
    Vm_r = optimizeCbModel(testModel,'min');
    FVA_r(i) = Vm_r.obj;
    
    fprintf('\b|\n');%for simple progress monitor
end
save('FVA.mat','FVA_f','FVA_r')
lowFluxRxns = model.rxns( FVA_f > 1e-9 & FVA_f < 1);
lowFluxRxns = union(lowFluxRxns, model.rxns( FVA_r < -1e-9 & FVA_r > -1));
save('lowFluxRxns.mat','lowFluxRxns')

