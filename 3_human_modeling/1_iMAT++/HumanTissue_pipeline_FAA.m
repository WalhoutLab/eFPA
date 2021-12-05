%% About
% this is just for technical accident related to cluster job scheduling.
% May be not used.

%% setup
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
%% Load model
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


% create missing field
model = creategrRulesField(model);
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
% Prepare flux thresholds (epsilon values)
load('epsilon_major10_side1_epsilon001_k01.mat'); % see walkthrough_generic.m for guidance on generating the epsilon values
if(any(epsilon_f<9e-6 | epsilon_r<9e-6))
    error('epsilon too small')
end

% [epsilon_f, epsilon_r,capacity_f,capacity_r] = makeEpsilonSeq(model, model.rxns, 0.1, 0.01)
model = changeRxnBounds(model,'EX_majorNutr',1000,'u');


% load the core set 
coreRxns = {'CYOOm3i','HMR_6911','HMR_6912','HMR_6914','HMR_6916','HMR_6918','HMR_6921'};

% all the libraries to model
libs = {'RNA'};
majoruptakes = setdiff(model.rxns(model.S(strcmp(model.mets,'majorNutr'),:)<0),{'EX_majorNutr'});
exrxns = setdiff(model.rxns(findExcRxns(model)),model.rxns(13102:end));

for zz = 1:length(libs)
    lib = libs{zz};
    wd = './input/humanModel/GeneCategories/';
    mkdir(['output/humanModel/',lib,'/OFD']);
    TPM = readtable('./input/humanModel/RNATissueMedian_log2TPM.xlsx');
    conditions = TPM.Properties.VariableNames(2:end);

    %% Step 2: FAA
    mkdir(['output/humanModel/',lib,'/FAA']);
    % Define the par cluster
    myCluster = parcluster('local');
    myCluster.NumWorkers = 128;
    saveProfile(myCluster);
    parpool(20,'SpmdEnabled',false);% adjust according to your computing environment

    for i = 14:16
	sampleName = conditions{i};
        load(['output/humanModel/',lib,'/OFD/',sampleName,'.mat']);
        % setup FAA inputs
        targetRxns = model.rxns;
        parforFlag = 1;
        % run FAA by calling:
        [levels_f,levels_r] = FAA_MILP(myCSM.MILP_PFD,model, myCSM.OFD, epsilon_f,epsilon_r,targetRxns, parforFlag);
        save(['output/humanModel/',lib,'/FAA/',conditions{i},'_levels.mat'],'levels_f','levels_r');
        fprintf('level table of %s is saved!\n',conditions{i});
    end

end
