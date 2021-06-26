function mytbl = listRxn(model,flux,myMet)
% This is a simple flux tracker to show the flux around a metabolite
% according to a flux distribution
%
% USAGE:
%
%    mytbl = listRxn(model,flux,myMet)
%
% INPUTS:
%    model:             input RECON2.2 model (COBRA model structure)
%    flux:              the flux distribution to inspect
%    myMet:             the metabolite of interest to list flux around
%
% OUTPUT:
%   mytbl:             	a table showing all the flux that produces or
%                       comsumes myMet
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020
myrxns = model.rxns(any(model.S(strcmp(model.mets,myMet),:),1));
myInds = any(model.S(strcmp(model.mets,myMet),:),1);
myfluxP = flux(ismember(model.rxns,myrxns)) .*  model.S(strcmp(model.mets,myMet),myInds)';
for i = 1:length(myrxns)
    mytbl(i,1) = myrxns(i);
    mytbl(i,2) = {flux(strcmp(model.rxns,myrxns{i}))};
    mytbl(i,3) = printRxnFormula_XL(model, myrxns{i},0);
    mytbl(i,4) = printRxnFormula(model, myrxns{i},0);
end
mytbl(:,5) = mat2cell(myfluxP,ones(length(myfluxP),1),1);
mytbl = sortrows(mytbl,5);
end