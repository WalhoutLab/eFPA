function mytbl = listFPtbl(model,solution)
model_irrev = convertToIrreversible(model);
% unify nomenclature - the nomenclature of reactions should be consistent 
% with the distance calculator
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); % reactions that have no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);% all "_f"
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
myrxns = solution.rxns;
mytbl(:,1) = myrxns;
mytbl(:,2) = printRxnFormula(model_irrev, myrxns,false);
mytbl(:,3) = num2cell(solution.full, 2);
mytbl(:,4) = num2cell(solution.weight,2);
mytbl(:,5) = num2cell(solution.weight .* solution.full, 2);
mytbl(:,6) = num2cell(solution.penalty, 2);
mytbl(:,7) = num2cell(solution.pDist', 2);
mytbl(:,8) = num2cell(solution.Dist', 2);
mytbl = sortrows(mytbl,5,'descend');
mytbl = [{'rxn','formula','flux','weight','allowance','penalty','pDistance','distance'};mytbl];
end