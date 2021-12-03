function model = loadYeatModel()
% setup the yeast model
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
model.geneNames(end+1) = {'VTC4'}; % gene name annotation for newly added gene
end
