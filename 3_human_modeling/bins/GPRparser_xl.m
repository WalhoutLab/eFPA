function parsedGPR = GPRparser_xl(model)
% Maps the GPR rules of the model to a specified format that is used by
% the model extraction methods 
%
% USAGE:
%   parsedGPR = GPRparser(model)
%
% INPUT:
%   model:       cobra model structure
%
% OUTPUT:
%   parsedGPR:   cell matrix containing parsed GPR rule
%
% AUTHORS: Thomas Pfau & Anne Richelle, May 2017
% Modified By: Xuhang Li, March 2020
%         - to fit into modified GPR parser "FormulaParser_xl"

parsedGPR = {};
fp = FormulaParser_xl();
for i = 1:numel(model.rxns)
    if ~isempty(model.rules{i})
        head = fp.parseFormula(model.rules{i});
        currentSets = head.getFunctionalGeneSets(model.genes)';
        parsedGPR{i}=currentSets;
    else
        parsedGPR{i}={''};
    end
end
parsedGPR=parsedGPR';
end