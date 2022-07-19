function HumanTissue_pipeline_fast(startTissueInd,endTissueInd)
startTissueInd = str2num(startTissueInd);
endTissueInd = str2num(endTissueInd);
%% Step 1: make the OFDs
% add paths
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

% correct a wrong GPR annotation that caused problem in the integration
model = changeGeneAssociation(model, 'HMR_4137',...
    'ENSG00000091140 and ENSG00000110435 and (ENSG00000131828 or ENSG00000163114) and ENSG00000150768 and ENSG00000168291');


% create missing field
model = creategrRulesField(model);
parsedGPR = GPRparser_xl(model);% Extracting GPR data from model
model.parsedGPR = parsedGPR;
% Prepare flux thresholds (epsilon values)
load('epsilon_major10_side1_epsilon001_k01.mat'); 
if(any(epsilon_f<9e-6 | epsilon_r<9e-6))
    error('epsilon too small. Check inputs!')
end

% ensure the exchange is unlimited
model = changeRxnBounds(model,'EX_majorNutr',1000,'u');


% load the core set 
coreRxns = {'CYOOm3i','HMR_6911','HMR_6912','HMR_6914','HMR_6916','HMR_6918','HMR_6921'};

% build the models 
libs = {'protein','RNA'};
majoruptakes = setdiff(model.rxns(model.S(strcmp(model.mets,'majorNutr'),:)<0),{'EX_majorNutr'});
exrxns = setdiff(model.rxns(findExcRxns(model)),model.rxns(13102:end));

for zz = 1:length(libs)
    lib = libs{zz};
    wd = './input/humanModel/GeneCategories/';
    mkdir(['output/humanModel/',lib,'/OFD']);
    TPM = readtable('./input/humanModel/RNATissueMedian_log2TPM.xlsx');
    conditions = TPM.Properties.VariableNames(2:end);

    % Set parameters
    modelType = 3; 
    speedMode = 2;
    %
    % First calculate and save the OFDs of all these conditions
    for i = startTissueInd:endTissueInd
        sampleName = conditions{i};
        load([wd,lib,'_refined_',sampleName,'.mat']);

        myCSM = struct(); % myCSM: my Context Specific Model
        [myCSM.OFD,...
        myCSM.PFD,...
        myCSM.N_highFit,...
        myCSM.N_zeroFit,...
        myCSM.minLow,...
        myCSM.minTotal,...
        myCSM.minTotal_OFD,...
        myCSM.MILP,...
        myCSM.MILP_PFD,...
        myCSM.HGenes,...
        myCSM.RLNames,...
        myCSM.OpenGene,...
        myCSM.latentRxn,...
        myCSM.Nfit_latent,...
        myCSM.wasteDW]...
        = IMATplusplus(model,epsilon_f,epsilon_r, ExpCateg, modelType,speedMode,...
        1e-5,true, true,0.05,1.1,[],[],[],1, coreRxns);
        save(['output/humanModel/',lib,'/OFD/',sampleName,'.mat'],'myCSM');
        eval([sampleName,' = myCSM;']);
        fprintf('OFD of %s is done!\n',sampleName);
        %%
        % inspect nutrient waste
        activeuptke = model.rxns(ismember(model.rxns,majoruptakes) & myCSM.OFD < -1e-7);
        wastedmass = 0;
        for j = 1:length(activeuptke)
            mymet = setdiff(model.mets((model.S(:,strcmp(model.rxns,activeuptke{j})) <0)),{'majorNutr'});
            myex = intersect(model.rxns(model.S(strcmp(model.mets,mymet),:)<0),exrxns);
            wastedmass = wastedmass + MolMass_ihuman(model.metFormulas{strcmp(model.mets,mymet)}) * myCSM.OFD(strcmp(model.rxns,myex{:}));
        end
        totalMass = sum(myCSM.OFD(ismember(model.rxns,{'EX_majorNutr'}))) * 10 * 100;
        wasteperc = wastedmass/totalMass;
        fprintf('wasted major percentage is %f!\n',wasteperc);
    end

    %% Step 2: FAA
    mkdir(['output/humanModel/',lib,'/FAA']);
    % Define the par cluster
    myCluster = parcluster('local');
    myCluster.NumWorkers = 128;
    saveProfile(myCluster);
    parpool(20,'SpmdEnabled',false);% adjust according to your computing environment

    for i = startTissueInd:endTissueInd
        sampleName = conditions{i};
        eval(['myCSM = ',sampleName,';']);
        % setup FAA inputs
        targetRxns = model.rxns;
        parforFlag = 1;
        % run FAA by calling:
        [levels_f,levels_r] = FAA_MILP(myCSM.MILP_PFD,model, myCSM.OFD, epsilon_f,epsilon_r,targetRxns, parforFlag);
        save(['output/humanModel/',lib,'/FAA/',conditions{i},'_levels.mat'],'levels_f','levels_r');
        fprintf('level table of %s is saved!\n',conditions{i});
    end
    delete(gcp('nocreate'))
end
end
