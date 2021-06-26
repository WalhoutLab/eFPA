%% custom reactions are already combined 
load('./../../../../yeastGEM.mat');
newyeast = model;
rxnTbl = readtable('rxn_yeast_full.tsv','FileType','delimitedtext','Delimiter','\t');
NaNrxns = rxnTbl(isnan(rxnTbl.StoiCoef),:); % they look protein names, we skip them
rxnTbl = rxnTbl(~isnan(rxnTbl.StoiCoef),:);

cmpTbl = readtable('comp_yeast.tsv','FileType','delimitedtext','Delimiter','\t');

revTbl = readtable('flux_dir_yeast_full.tsv','FileType','delimitedtext','Delimiter','\t');

grTbl = readtable('rxn_par_yeast_full.tsv','FileType','delimitedtext','Delimiter','\t');

metTbl = readtable('spec_yeast_full.tsv','FileType','delimitedtext','Delimiter','\t');

model = createModel();
rxns = unique(rxnTbl.ReactionID);
for i = 1:length(rxns)
    metlist = rxnTbl.Metabolite(ismember(rxnTbl.ReactionID,rxns(i)));
    [A B] = ismember(metlist,metTbl.SpeciesID);
    complist = metTbl.Compartment(B(A));
    for j = 1:length(metlist)
        metlist{j} = [metlist{j},'[',complist{j},']'];
    end
    rev = revTbl.Reversible(ismember(revTbl.ReactionID,rxns(i)));
    rev = strcmp(rev{:},'true');
    % according to simmer code: reversibleRx$reversible = ifelse(reversibleRx$modelBound == "greaterEqual", 1, 0) # directionality is coded as -1: irreversibly backward, 0: reversible, 1: irreversibly forward
    if strcmp(revTbl.FluxBound{ismember(revTbl.ReactionID,rxns(i))},'greaterEqual')
        rev = false;
    end
    if ~rev
        lb = 0;
    else
        lb = -1000;
    end
    
    grrule = grTbl.Enzymes{strcmp(grTbl.ReactionID,rxns{i})};
    stolist = rxnTbl.StoiCoef(ismember(rxnTbl.ReactionID,rxns(i)));
    model = addReaction(model,rxns{i},'metaboliteList',metlist,'stoichCoeffList',stolist, 'reversible',rev,...
        'lowerBound',lb,'upperBound',1000,'geneRule',grrule,'printLevel',0);
end
%% additional reversibility from themo
model.lb(ismember(model.rxns,{'r_0659','r_0470'})) = 0;
%% add met name
metid = model.mets;
metid = regexprep(metid,'\[c_..\]$','');
[A B] = ismember(metid, metTbl.SpeciesID);
model.metNames = metTbl.SpeciesName(B(A));
%% set up boundarys and exchanges 
% first get the default bounds from newest yeast 
openEx = newyeast.rxns(findExcRxns(newyeast) & newyeast.lb < 0);
[A B] = ismember(openEx,newyeast.rxns);
openBounds = newyeast.lb(B(A));
% second convert the bound met to EX rxns and set the bounds
bondMet = model.mets(cellfun(@(x) ~isempty(regexp(x,'\[c_01\]$', 'once')),model.mets));
bondRxn = model.rxns(any(model.S(ismember(model.mets,bondMet),:),1));
for i = 1:length(bondRxn)
    model.S(ismember(model.mets,bondMet),strcmp(model.rxns,bondRxn{i})) = 0;
    if any(strcmp(openEx,bondRxn{i}))
        model.lb(strcmp(model.rxns,bondRxn{i})) = openBounds(strcmp(openEx,bondRxn{i}));
    else
        model.lb(strcmp(model.rxns,bondRxn{i})) = 0;
    end
end
%% set up others
% set ATP maintenance
% add maintance rxns
model = addReaction(model,['r_4046'],'reactionFormula','s_0434[c_03] + s_0803[c_03]  -> s_0394[c_03] + s_0796[c_06] + s_1322[c_03]','geneRule', '','printLevel',1);
model = changeRxnBounds(model,'r_4046',1,'b'); % maintance 

model = buildRxnGeneMat(model);
model = creategrRulesField(model);
%% test
model = changeObjective(model,'r_2111');
optimizeCbModel(model)
save('./../../../../simmer_model.mat','model');
%% check consistency
model_ori = readSBML('./../input/YeastJoshua/originalDataTbl/yeast_6.00_model.xml',1000);
rxn2check = readtable('./../input/YeastJoshua/originalDataTbl/rxnID2check.xlsx');
%%
for i = 2: length(rxn2check.rxns)
    fm1 = printRxnFormula_XL(newyeast,rxn2check.rxns{i},false);
    fm2 = printRxnFormula_XL(model,rxn2check.rxns{i},false);
    if ~strcmp(fm1{:},fm2{:})
        rxn2check.rxns{i}
        fm1{:}
        fm2{:}
        break;
    end
end
%%
for i = 219: length(rxn2check.rxns)
    fm1 = printGPRForRxns(newyeast,rxn2check.rxns{i},false);
    fm2 = printGPRForRxns(model,rxn2check.rxns{i},false);
    if ~strcmp(fm1{:},fm2{:})
        rxn2check.rxns{i}
        fm1{:}
        fm2{:}
        break;
    end
end
% broken indexes are: 
% r_0005 r_0148 r_0226 r_0250 r_0364 r_0366 r_0534 r_0542 r_0886 r_0887
% r_0888 r_0893 r_0916 r_1021 r_1084 r_1166 r_1172 r_1574 r_1667 r_2034
% r_3526
% 10% GPR difference 
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

    
    
    
    
    
    
    