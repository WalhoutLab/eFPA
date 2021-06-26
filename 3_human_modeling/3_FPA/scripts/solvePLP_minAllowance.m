function [obj, solution] = solvePLP_minAllowance(model,w, w_label, obj_rxn,default_w, minAllowance,minMax) 
% form and solve a special Linear Problem (LP) we named "Penalized Linear
% Problem (PLP)". It creates a LP for the following:
% maximize TargetRxn
%       s.t. (1) w * v <= a (by default a = 1)
%            (2) other LP constriants
%       w: weight vector
%       v: flux vector
%       a: flux allowance
%
% USAGE:
%
%    obj = solvePLP(model,w, w_label, obj_rxn,default_w, minMax) 
%
% INPUTS:
%    model:             input model (COBRA model structure, irreversible model)
%    w:                 the weight vector. Not required to be equal to number of reactions in the model (reactions that don't
%                       have a specified weight will be of the default
%                       weight)
%    w_label:           reaction ID of the weight vector. Should be in
%                       equal length of w
%    obj_rxn:           the objective reaction to maximize
%    default_w:         the default weight if a reaction is not specified
%                       for a weight in w
%    minMax:            doing maximization or minimization (default:
%                       maximization)
%
%
% OUTPUT:
%   obj:                the objective value of PLP, aka, flux potential.
% OPTIONAL OUTPUTS:
%   solution:           the FPA solution
%
% `Yilmaz et al. (2020). Final Tittle and journal.
%
% .. Author: - (COBRA implementation) Xuhang Li, Mar 2020

if (nargin < 7)
    minMax = 'max';
end
%% step1: match with the weight(w) input
% normally the weight should already be aligned with the rxns in the model,
% but solvePLP function also allows input partial weight list with a defined
% label set.So, we first need to match the weight to reactions.
% match the weight
w_full = default_w .* ones(length(model.rxns),1);
[A B] = ismember(model.rxns, w_label);
w_full(A) = w(B(A));
%% add the penalty constriants 
model.S(end+1,:) = w_full';
model.csense(end+1) = 'L';
model.b(end+1) = minAllowance; %flux allowance.
% also update other related fields (could be skiped)
if isfield(model,'BiGG')
    model.BiGG(end+1) = {'NA'};
end
if isfield(model,'metCharges')
    model.metCharges(end+1) = 0;
end
if isfield(model,'metFormulas')
    model.metFormulas(end+1) = {'NA'};
end
if isfield(model,'metNames')
    model.metNames(end+1) = {'Nonsense Metabolites which is a constraint'};
end
if isfield(model,'mets')
    model.mets(end+1) = {['NonMetConst',num2str(length(model.mets))]};
end
% solve the LP
model = changeObjective(model,obj_rxn);
solution = optimizeCbModel(model,minMax);
if solution.stat == 1
    obj = solution.obj;
else
    obj = NaN;
end
solution.weight = w_full;%we save the weight vector in solution output for tracking purposes
