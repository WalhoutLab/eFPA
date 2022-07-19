%% About 
% this is to write out the model (to plain text tables) for the distance
% calculator program 

%% Add required path 
addpath ./input/
addpath ./scripts/
addpath ./../input/
addpath ./../bins/
initCobraToolbox(false);
%% 1. Prepare the model
load('input/ihuman_COBRA.mat');
mkdir('./../distance_inputs');
%% 2. write
writematrix(model.S,'distance_inputs/Smatrix_regular.txt');
writecell(model.rxns,'distance_inputs/reactions_regular.txt');
writecell(model.mets,'distance_inputs/metabolites_regular.txt');
writematrix(model.lb,'distance_inputs/LB_regular.txt');
writematrix(model.ub,'distance_inputs/UB_regular.txt');
byProducts = {'CO2';'AMP';'NADP+';'NADPH';'PPi';'O2';'NADH';'NAD+';
            'Pi';'ADP';'CoA';'ATP';'H2O';'H+';'GTP';'GDP';
            'Electron Transfer Flavoprotein Reduced';'Electron Transfer Flavoprotein Oxidized';
            'L-carnitine';'FAD';'FADH2';'Na+'};% the typical set of byproducts
% Add compartment label to byproducts
byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,' \[.+\]$',''),model.metNames, 'UniformOutput',false),byProducts));
writecell(byProducts,'distance_inputs/byproducts_regular.txt');
