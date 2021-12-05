%% About
% making the epsilons for iMAT++ integration
%%
addpath ~/cobratoolbox/% your cobra toolbox path
addpath /share/pkg/gurobi/900/linux64/matlab/% the gurobi path
addpath ./../bins/
addpath ./../input/
addpath scripts/
initCobraToolbox(false);
% Load model
load('input/ihuman_serum.mat');
model.subSystems = [model.subSystems{:}]';
% the default constraints are unlimited, so we dont change anything

% ```special treatments for iHuman1```
% Remove parentathsis in the reaction ID (which will be changed in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% change the name starting with numbers 
model.rxns = regexprep(model.rxns,'(^[0-9])','x$1');

model = changeRxnBounds(model,'EX_majorNutr',10,'u');
model = changeRxnBounds(model,'EX_sideNutr',1,'u');
reactions = model.rxns;
epsilon0 = 0.01;
k = 0.1;
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

myCluster = parcluster('local');
myCluster.NumWorkers = 128;
saveProfile(myCluster);
parpool(20,'SpmdEnabled',false); % we run this code on a 40-core lab server

epsilon_f = epsilon0 * ones(length(reactions),1);
epsilon_r = epsilon0 * ones(length(reactions),1);
capacity_f = zeros(length(reactions),1);
capacity_r = zeros(length(reactions),1);
FVA_f = zeros(length(reactions),1);
FVA_r = zeros(length(reactions),1);
environment = getEnvironment();
tol = 1e-8; % the default flux numerical tolerance for capacity output
fprintf(['\n' repmat('.',1,length(reactions)) '\n\n']);
parfor i = 1: length(reactions)
    restoreEnvironment(environment);
    testModel = model;
    testModel = changeObjective(testModel,reactions(i));
    Vm_f = optimizeCbModel(testModel, 'max');
    FVA_f(i) = Vm_f.obj;
    if Vm_f.obj*k < epsilon0 && Vm_f.obj*k > 1e-10 %avoid numeric error
        epsilon_f(i) = Vm_f.obj * k;
    end
    if Vm_f.obj > tol
        capacity_f(i) = true;
    else
        capacity_f(i) = false;
    end

    Vm_r = optimizeCbModel(testModel,'min');
    FVA_r(i) = Vm_r.obj;
    if -Vm_r.obj*k < epsilon0 && Vm_r.obj*k < -1e-10
        epsilon_r(i) = -Vm_r.obj * k;
    end
    if Vm_r.obj < -tol
        capacity_r(i) = true;
    else
        capacity_r(i) = false;
    end
    fprintf('\b|\n');%for simple progress monitor
end
save('epsilon_major10_side1_epsilon001_k01.mat','epsilon_f','epsilon_r','capacity_f','capacity_r','FVA_f','FVA_r')
