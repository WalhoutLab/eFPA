function mytbl = listFPtbl(model,solution)
model_irrev = convertToIrreversible(model);
% unify nomenclature - the nomenclature of reactions should be consistent 
% with the distance calculator
% tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); % reactions that have no suffix
% model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);% all "_f"
% model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
myrxns = solution.rxns;
mytbl(:,1) = myrxns;
mytbl(:,2) = model_irrev.rxnNames;
[A B] = ismember(myrxns, model_irrev.rxns);
mytbl(:,3) = model_irrev.subSystems(B(A));
mytbl(:,4) = printRxnFormula_XL(model_irrev, myrxns,false);
mytbl(:,5) = num2cell(solution.full, 2);
mytbl(:,6) = num2cell(solution.weight,2);
mytbl(:,7) = num2cell(solution.weight .* solution.full, 2);
mytbl(:,8) = num2cell(solution.penalty, 2);
mytbl(:,9) = num2cell(solution.pDist', 2);
mytbl(:,10) = num2cell(solution.Dist', 2);
mytbl(:,11) = printGPRForRxns(model_irrev, myrxns,false);

mytbl = sortrows(mytbl,7,'descend');
mytbl = [{'rxn','rxnName','subsys','formula','flux','weight','allowance','penalty','pDistance','distance','GPR'};mytbl];
end