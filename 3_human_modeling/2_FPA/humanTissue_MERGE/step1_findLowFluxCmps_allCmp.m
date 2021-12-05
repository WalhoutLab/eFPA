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

load('allCmp_iHumanName.mat');
cmpOrder = {'[Cytosol]','[Mitochondria]','[Inner mitochondria]',...
            '[Peroxisome]','[Lysosome]','[Golgi apparatus]',...
            '[Endoplasmic reticulum]','[Nucleus]','[Extracellular]'}; % some met like sucrose and pectin only exist in ex cellular

% create missing field
model = creategrRulesField(model);
cmps = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
cmps = setdiff(cmps,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
environment = getEnvironment();
fprintf(['\n' repmat('.',1,length(cmps)) '\n\n']);
parfor i = 1: length(cmps)
    restoreEnvironment(environment);
    
    % determine the avaiable compartment 
    myMet = cmps{i};
    metInd = cellfun(@(x) ~isempty(regexp(x,['^',regexptranslate('escape',myMet),' \[(\w|\s)*\]$'],'once')),model.metNames);
    allCmp = regexp(model.metNames(metInd),'\[(\w|\s)*\]$','match');
    allCmp = unique([allCmp{:}]);
    targetCmp = cmpOrder(ismember(cmpOrder,allCmp));
    targetCmp = targetCmp{1};
    metID = model.mets(strcmp(model.metNames,[myMet,' ',targetCmp]));
    % add demand rxn
    testModel = addReaction(model,'targetDMN','metaboliteList',metID,'stoichCoeffList',-1,'geneRule', 'NA','printLevel',0);  
    testModel = changeObjective(testModel,'targetDMN');
    Vm_f = optimizeCbModel(testModel, 'max');
    FVA_f(i) = Vm_f.obj;
    
    Vm_r = optimizeCbModel(testModel,'min');
    FVA_r(i) = Vm_r.obj;
    
    fprintf('\b|\n');%for simple progress monitor
end
save('FVA_cmp_all.mat','FVA_f','FVA_r')
lowFluxCmps = cmps( FVA_f > 1e-9 & FVA_f < 1);
lowFluxCmps = union(lowFluxCmps, cmps( FVA_r < -1e-9 & FVA_r > -1));
save('lowFluxCmp_all.mat','lowFluxCmps')

