function [expressionRxns] = mapExpressionToReactions_xl(model, expressionData)                                          
% Determines the expression data associated to each reaction present in
% the model 
%
% USAGE:
%    [expressionRxns parsedGPR] = mapExpressionToReactions_xl(model, expressionData) 
%
% INPUTS:
%	model                   model strusture
%	expressionData          mRNA expression data structure
%       .gene               	cell array containing GeneIDs in the same
%                               format as model.genes
%       .value                  value should be qudra-level for the
%                               integration. high exp is 3; moderate exp is 2; low is 1; zero
%                               is 0
% OUTPUTS:
%   expressionRxns:         reaction expression, corresponding to model.rxns.
%   parsedGPR:              cell matrix containing parsed GPR rule
%
% .. Authors:
%       - Anne Richelle, May 2017 - integration of new extraction methods 
% Modified By: Xuhang Li
%           - slightly modified for saving time: instead of parsing GPR,
%               use the preparsed GPR in the model

parsedGPR = model.parsedGPR;% Extracting GPR data from model
% Find wich genes in expression data are used in the model
[gene_id, gene_expr] = findUsedGenesLevels(model,expressionData);
% Link the gene to the model reactions
expressionRxns = selectGeneFromGPR_xl(model, gene_id, gene_expr, parsedGPR);
