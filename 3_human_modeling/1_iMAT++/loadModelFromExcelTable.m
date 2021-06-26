%% to make the corbra model from an excel table
%% the generic C. elegans model
addpath ./../bins/
% initCobraToolbox;
% the <- in rxns table should be change to <==> manually
% ub/lb header to "Lower bound	Upper bound"
% change met header
% add two psudo met to the met table
% the rxns table and met table should be merged into one excel table
model = iCEL2model_xl('./originalModel.xlsx');
% test model
model = changeObjective(model,'BIO0010');
opt = optimizeCbModel(model,'max')
model = changeObjective(model,'RCC0005');
optimizeCbModel(model,'max')
% save model
model = changeObjective(model,'BIO0010');
model.description = 'iCEL1314';
save('iCEL1314.mat','model');

%% make the dual model with dynamic controls
model = iCEL2model_xl('./dualModel.xlsx');
sideR=0.02;
sideR_ind=0.005;
storageR=0.01;
bacMW=966.28583751;
% dividedBy=100; this is in the stoicheomitry matrix
% add individual met control and scale by the factor
sideRxns =  model.rxns(model.S(strcmp(model.mets,'sideMet[e]'),:)~=0);
sideRxns(strcmp(sideRxns,'EXC9998_E')) = [];
for i = 1:length(sideRxns)
    model.S(end+1,:) = zeros(1,length(model.rxns));
    model.S(end, strcmp(sideRxns{i},model.rxns)) = model.S(strcmp(model.mets,'sideMet[e]'),strcmp(sideRxns{i},model.rxns)); 
    model.S(end, strcmp('EXC0050_L',model.rxns)) = sideR_ind*bacMW*0.01;%the 0.01 here is the 1/devidedBy 
    model.csense(end+1) = 'L';
    model.b(end+1) = 0;
    model.BiGG(end+1) = {'NA'};
    model.metCharges(end+1) = 0;
    model.metFormulas(end+1) = {'NA'};
    model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
    model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
end

% add overall side/storage control
model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9999_E',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050_L',model.rxns)) = storageR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model.S(end+1,:) = zeros(1,length(model.rxns));
model.S(end, strcmp('EXC9998_E',model.rxns)) = 1; 
model.S(end, strcmp('EXC0050_L',model.rxns)) = sideR*bacMW*0.01; 
model.csense(end+1) = 'L';
model.b(end+1) = 0;
model.BiGG(end+1) = {'NA'};
model.metCharges(end+1) = 0;
model.metFormulas(end+1) = {'NA'};
model.metNames(end+1) = {'Pseudo-metabolite that represents a constraint'};
model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};

model = changeObjective(model,'BIO0010_I');
opt = optimizeCbModel(model,'max')
model = changeObjective(model,'RCC0005_X');
optimizeCbModel(model,'max')
% save the model
model = changeObjective(model,'BIO0010_I');
model.description = 'iCEL1314 - dual model';
save('Tissue.mat','model');
%% model test
ref = readtable('dualModelTable.tsv','FileType','text');
%%
rxn = 'BIO0010_I';
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);
%test = changeRxnBounds(model,'DMN0033',1,'u');
model = changeObjective(model,rxn);
opt = optimizeCbModel(model,'max');
test = changeRxnBounds(model,rxn,opt.obj-1e-5,'l');
minflux = minimizeModelFlux_XL(test);
sum(abs(ref.(rxn))) - sum(abs(minflux))
