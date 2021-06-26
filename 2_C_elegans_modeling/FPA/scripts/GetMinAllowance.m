function [allowance, Vmax, solution] = GetMinAllowance(model,w, w_label, obj_rxn,default_w, minMax) 
% form and solve a special Linear Problem (LP) to calculate the minimum
% flux allowance for a target reaction. It creates a LP for the following:
% minimize a
%       s.t. (1) w * v <= a
%            (2) v(obj_rxn) >= v.max(obj_reaction)
%            (2) other LP constriants
%       w: weight vector
%       v: flux vector
%       a: flux allowance
%
% USAGE:
%
%    obj = GetMinAllowance(model,w, w_label, obj_rxn,default_w, minMax) 
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
%   obj:                the minimum flux allowance for the objective
%                       reaction
%
% `Yilmaz et al. (2020). Final Tittle and journal.
%
% .. Author: - (COBRA implementation) Xuhang Li, Mar 2020

if (nargin < 6)
    minMax = 'max';
end
if (nargin < 5)
    default_w = 0; %the default weight for a reaction if it is not included in any weight input; 0 means not contributing to the total flux allowance
end

%% step1: prepare the FBA objective value
% first perform FBA to calculate the maximum flux, Vm
model = changeObjective(model,obj_rxn);
solVm = optimizeCbModel(model,minMax);
if strcmp(minMax, 'max')
    model.lb(strcmp(model.rxns,obj_rxn)) = solVm.obj-1e-9;
else
    model.ub(strcmp(model.rxns,obj_rxn)) = solVm.obj+1e-9;
end
%% step2: calculate the minimum flux allowance
% match with the weight(w) input
% normally the weight should already be aligned with the rxns in the model,
% but solvePLP function also allows input partial weight list with a defined
% label set.So, we first need to match the weight to reactions.
% match the weight
w_full = default_w .* ones(length(model.rxns),1);
[A B] = ismember(model.rxns, w_label);
w_full(A) = w(B(A));
%% create the allowance objective function
model.c = w_full;
% solve the LP 
solution = optimizeCbModel(model,'min');
if solution.stat == 1
    allowance = solution.obj;
    Vmax = solution.full(strcmp(model.rxns,obj_rxn));
else
    allowance = NaN;
    Vmax = NaN;
    solution = {''};
    warning('The LP is infeasible. Please check what is going wrong!');
    fprintf('The infeasible RXN is %s\n',obj_rxn);
end
