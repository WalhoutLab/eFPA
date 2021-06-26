function mytbl = listFPtbl(model,solution)
model_irr = convertToIrreversible(model);
myrxns = model_irr.rxns;
mytbl(:,1) = myrxns;
mytbl(:,2) = printRxnFormula_XL(model_irr, myrxns,false);
mytbl(:,3) = num2cell(solution.full, 2);
mytbl(:,4) = num2cell(solution.weight,2);
mytbl(:,5) = num2cell(solution.weight .* solution.full, 2);
mytbl(:,6) = num2cell(solution.penalty, 2);
mytbl(:,7) = num2cell(solution.pDist', 2);
mytbl(:,8) = num2cell(solution.Dist', 2);
mytbl = sortrows(mytbl,5,'descend');
mytbl = [{'rxn','formula','flux','weight','allowance','penalty','pDistance','distance'};mytbl];
end