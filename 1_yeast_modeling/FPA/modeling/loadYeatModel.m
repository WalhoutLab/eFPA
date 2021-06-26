function model = loadYeatModel()
load('./../../input/GEMs/yeastGEM.mat');
% add polyphosphate kinase
model = addReaction(model,['cstm_PPK'],'reactionFormula','s_1322[c] + s_0434[c] <==> s_0394[c] + s_4018[c]','geneRule', 'YJL012C','printLevel',1,'reversible',false,'lowerBound',0);
% prohibit the spontenous polymerization
model.lb(strcmp(model.rxns,'r_4333')) = 0;

% add orotate secretion
model = addReaction(model,['cstm_s_1269_TCE'],'reactionFormula','s_1269[c] <==> s_1269[e]','geneRule', '','printLevel',1);
model = addReaction(model,['cstm_s_1269_EX'],'reactionFormula','s_1269[e] <==> ','geneRule', '','printLevel',1,'lowerBound',0);

% set ATP maintenance and misc
model = changeRxnBounds(model,'r_4046',1,'b'); % maintance 
model = addReaction(model,['cstm_misc_ATPhy'],'reactionFormula','s_0434[c] + s_0803[c]  -> s_0394[c] + s_0794[c] + s_1322[c]','geneRule', '','printLevel',1);

model = buildRxnGeneMat(model);
model = creategrRulesField(model);
model.geneNames(end+1) = {'VTC4'};
%% check consistency
% model_ori = readSBML('./../input/YeastJoshua/originalDataTbl/yeast_6.00_model.xml',1000);
% rxn2check = readtable('./../input/YeastJoshua/originalDataTbl/rxnID2check.xlsx');
% %%
% for i = 233: length(rxn2check.rxns)
%     fm1 = printRxnFormula_XL(model_ori,rxn2check.rxns{i},false);
%     fm2 = printRxnFormula_XL(model,rxn2check.rxns{i},false);
%     if ~strcmp(fm1{:},fm2{:})
%         rxn2check.rxns{i}
%         fm1{:}
%         fm2{:}
%         break;
%     end
% end
% broken indexes are:
%   old==>new
% 'r_1099' ==> 'r_2132'
% 'r_4040' ==> 'cstm_PPK'
% 'r_4042' ==> 'cstm_misc_ATPhy'
% 'r_4043' ==> 'cstm_s_1269_TCE'
%% media composition
% Total nutrient = {'r_2005','r_1714','r_1654','r_2060','r_2090','r_1899'};
% freeEx = {'r_1672','r_2100','r_1832','r_1992'};
%% the following nutrients need to be set manually
% % phosphate exchange
% model.lb(strcmp(model.rxns,'r_2005')) = r_1244;
% % glucose exchange
% model.lb(strcmp(model.rxns,'r_1714')) = r_1166;
% % ammonium exchange 
% model.lb(strcmp(model.rxns,'r_1654')) = r_1115;
% % uracil 
% model.lb(strcmp(model.rxns,'r_2090')) = r_1272;
% % leucine
% model.lb(strcmp(model.rxns,'r_1899')) = r_1211;
%% write the distance input
% writematrix(model.S,'./../distance_inputs/Smatrix_regular.txt');
% writecell(model.rxns,'./../distance_inputs/reactions_regular.txt');
% writecell(model.mets,'./../distance_inputs/metabolites_regular.txt');
% writematrix(model.lb,'./../distance_inputs/LB_regular.txt');
% writematrix(model.ub,'./../distance_inputs/UB_regular.txt');
% byProducts = {'carbon dioxide';'AMP';'NADP(+)';'NADPH';'diphosphate';'oxygen';'NADH';'NAD';'phosphate';'ADP';'coenzyme A';'ATP';'H2O';'H+';'GTP';'GDP';'(R)-carnitine';'FAD';'FADH2'};
% % add compartment label to byproducts
% byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,'\s\[.*\]$',''),model.metNames, 'UniformOutput',false),byProducts));
% writecell(byProducts,'./../distance_inputs/byproducts_regular.txt');
end
