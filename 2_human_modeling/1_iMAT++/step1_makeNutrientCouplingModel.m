%% About
% reconstruct the modified human model that has the nutrient coupling
% system (see methods for details)
%% to make the corbra model from an excel table
% Load model
load('input/ihuman_COBRA.mat');
model.subSystems = [model.subSystems{:}]';
%% adjust the S matrix scale
maxCoef = 100;
[x y] = find(abs(model.S)>maxCoef);
rxns2adj = unique(model.rxns(y));
for i = 1:length(rxns2adj)
    maxS = max(abs(model.S(:,strcmp(model.rxns,rxns2adj{i}))));
    scaleFactor = maxCoef / maxS;
    model.S(:,strcmp(model.rxns,rxns2adj{i})) = model.S(:,strcmp(model.rxns,rxns2adj{i})) .* scaleFactor;
end
% display the maximum stoichemitry in the adjusted matrix
max(max(abs(model.S)))
%% add side and major coupling constraints
fluxBurdenCoef_major = 100; % one unit of uptake flux will always import 100ug of mass 
fluxBurdenCoef_side = 100; % one unit of uptake flux will always import 100ug of mass 

largeFluxOffsetFactor = 10;
minCoef = 1e-5;

% fix some missing formula 
model.metFormulas(strcmp(model,'hdl_hs_s[s]')) = {'C4393H7057N1219O1349S22P1R5'};
model.metFormulas(strcmp(model,'idl_hs_s[s]')) = {'C12507.5H19765.5N3294.5O3741.5S57R12'};
model.metFormulas(strcmp(model,'ldl_hs_s[s]')) = {'C46525H73399N12226O13915S208P2R7'};
model.metFormulas(strcmp(model,'vldl_hs_s[s]')) = {'C5023.2H7917N1296.4O1523.2S22.4P2R19'};
model.metFormulas(strcmp(model,'chylo_hs_s[s]')) = {'C27538H43557N7330O8284S126R3'};
% some reaction has no MW and cannot be taken up, we skip these to avoid
% error
skipRxns = {'HMR_10029','HMR_10030','HMR_10031','HMR_10032'};

% freeRxns = {'EX_HC02154[e]','EX_HC02199[e]','EX_xolest_hs[e]','HMR_10024','HMR_10029',...
%               'HMR_9700','HMR_9687','HMR_10030','HMR_10031','HMR_10032','HMR_9162','HMR_9202','HMR_9204','HMR_9686'}; 
% these met cannot be taken up and do not have a defined formula, we leave
% them for free exchange ==> we decide to put artificial MW (protein
% average) to leave then in minor nutrient set, so that they wouldnt get too
% large flux to cause weird effect (like generate NADPH)


majorMets = readtable('input/humanModel/MajorBloodNutrients.xlsx');
inogMets = readtable('input/humanModel/MajorBloodNutrients.xlsx','Sheet','inorganicMet');
% first, block all uptake direction but leave free the inorganic ones
allEx = model.rxns(findExcRxns(model));
model = changeRxnBounds(model,allEx,0,'l');
model = changeRxnBounds(model,inogMets.RxnID,-1000,'l');
model = changeRxnBounds(model,skipRxns,-1000,'l');
% then, create the uptake reactions for all organic contents 
ToUptakeMajor = majorMets.RxnID;
ToUptakeMinor = setdiff(allEx,union(majorMets.RxnID,inogMets.RxnID));
ToUptakeMinor = setdiff(ToUptakeMinor,skipRxns);
% add the two proxy metabolites
model = addMetabolite(model, 'majorNutr', 'Major Nutrient in Blood Serum');
model = addMetabolite(model, 'sideNutr', 'Side Nutrient in Blood Serum');
% add the uptake reactions
for i = 1:length(ToUptakeMajor)
    myMet = model.mets{model.S(:,strcmp(model.rxns,ToUptakeMajor{i}))<0};
    WM = MolMass_ihuman(model.metFormulas{strcmp(model.mets,myMet)});
    coef = fluxBurdenCoef_major ./ WM;
    if coef < minCoef
        scaleFactor = minCoef / coef;
    else
        scaleFactor = 1;
    end
    model = addReaction(model,['UPK',myMet],'metaboliteList',{myMet,'majorNutr'},'stoichCoeffList',[-coef*scaleFactor, -1*scaleFactor],'geneRule', 'NA','lowerBound',-1000,'upperBound',0,'printLevel',0);
end
% each unit of master nutrient exchange will secrect 10 majorNutr (or
% sideNutr), making 1g of weight. This is to avoid the flux of these two
% master exchange becomes too large compared with the flux scale of others
model = addReaction(model,'EX_majorNutr','metaboliteList',{'majorNutr'},'stoichCoeffList', -largeFluxOffsetFactor,'geneRule', 'NA','lowerBound',0,'upperBound',1000,'printLevel',0);

for i = 1:length(ToUptakeMinor)
    myMet = model.mets{model.S(:,strcmp(model.rxns,ToUptakeMinor{i}))<0};
    WM = MolMass_ihuman(model.metFormulas{strcmp(model.mets,myMet)});
    coef = fluxBurdenCoef_side ./ WM;
    if coef < minCoef
        scaleFactor = minCoef / coef;
    else
        scaleFactor = 1;
    end
    model = addReaction(model,['UPK',myMet],'metaboliteList',{myMet,'sideNutr'},'stoichCoeffList',[-coef*scaleFactor, -1*scaleFactor],'geneRule', 'NA','lowerBound',-1000,'upperBound',0,'printLevel',0);
end
model = addReaction(model,'EX_sideNutr','metaboliteList',{'sideNutr'},'stoichCoeffList', -largeFluxOffsetFactor,'geneRule', 'NA','lowerBound',0,'upperBound',1000,'printLevel',0);

max(max(abs(model.S)))
a = abs(model.S);
min(min(a(a>0)))
%% save model for FVA
save('input/ihuman_serum_small_flux_unadjusted.mat','model');
%% run the following part after running FVA on the unadjusted model;
%% adjust the stoichemitry of very small flux reactions 
load('epsilon_FVA.mat');
rxns2adj = model.rxns(FVA_f > 0 & FVA_f < 1e-4);
rxns2adj = union(rxns2adj,model.rxns(FVA_r < 0 & FVA_r > -1e-4));
for i = 1:length(rxns2adj)
    minFluxf = FVA_f(strcmp(model.rxns,rxns2adj{i}));
    minFluxr = FVA_r(strcmp(model.rxns,rxns2adj{i}));
    minFlux = [];
    if minFluxf > 0 && minFluxf < 1e-4
        minFlux = [minFlux,minFluxf];
    end
    if minFluxr < 0 && minFluxr > -1e-4
        minFlux = [minFlux,minFluxr];
    end
    minFlux = min(abs(minFlux));
    adjFactor = minFlux /1e-4;
    model.S(:,strcmp(model.rxns,rxns2adj{i})) = model.S(:,strcmp(model.rxns,rxns2adj{i})) .* adjFactor;
end
%% TEST AND SAVE MODEL
model = changeObjective(model,'biomass_human');
opt = optimizeCbModel(model,'max')
model = changeObjective(model,'HMR_6916');
optimizeCbModel(model,'max')
% save the model
model = changeObjective(model,'biomass_human');
model.description = 'iHuman - serum';
save('input/ihuman_serum.mat','model');
%%
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
test = changeObjective(model,'HMR_6693');
test = changeRxnBounds(test,'EX_majorNutr',10,'u');
test = changeRxnBounds(test,'EX_sideNutr',1,'u');
optimizeCbModel(test)