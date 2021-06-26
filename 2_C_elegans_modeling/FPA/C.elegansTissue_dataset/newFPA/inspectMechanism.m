initCobraToolbox
expressionTbl = readtable('expressionTable.tsv','FileType','text','ReadRowNames',true);
tissueLabel = expressionTbl.Properties.VariableNames;
parpool(2)
load('Tissue.mat');

%%
targetObj = {'TCE0453'};
targetDir = 'f';
%%
[relFP_f,relFP_r, FluxPotential_solutions_X,FluxPotential_solutions_I] = run_a_case_originalMERGE(targetObj,1);
%%
[relFP_f,relFP_r, FluxPotential_solutions_X,FluxPotential_solutions_I] = run_a_case_wtdDist(targetObj,1);
%%
[relFP_f,relFP_r, FluxPotential_solutions_X,FluxPotential_solutions_I] = run_a_case_wtdDist_exp_decay(targetObj,1);
%%
figure(1)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,eval(['relFP_',targetDir]));
%% {'Gonad'}    {'Glia'}    {'Intestine'}    {'Pharynx'}    {'Hypodermis'}    {'Neurons'}    {'Body_wall_muscle'}
qryTissue = 'Glia';

if ~strcmp(qryTissue, 'Intestine')
    mytbl = listFPtbl(model,FluxPotential_solutions_X{find(strcmp(tissueLabel,qryTissue))}{strcmp(targetDir,'r')+1});
else
    mytbl = listFPtbl(model,FluxPotential_solutions_I{find(strcmp(tissueLabel,qryTissue))}{strcmp(targetDir,'r')+1});
end

%% titrate the distance order as needed
[relFP_f,relFP_r, FluxPotential_solutions_X,FluxPotential_solutions_I] = Titrate_A_Met_originalMERGE(targetObj,0:0.5:10);
%%
[relFP_f,relFP_r, FluxPotential_solutions_X,FluxPotential_solutions_I] = Titrate_A_Met_wtdDist(targetObj,0:0.5:10);
%%
[relFP_f,relFP_r, FluxPotential_solutions_X,FluxPotential_solutions_I] = Titrate_A_Met_wtdDist_exp_decay(targetObj,[0:0.5:10,15, 20, 50]);