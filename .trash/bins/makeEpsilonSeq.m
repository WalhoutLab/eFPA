function [epsilon_f, epsilon_r,capacity_f,capacity_r] = makeEpsilonSeq(model, reactions, epsilon0, k)
% Generate the epsilon sequence for performing IMAT or IMAT++ integration.
% The epsilons are determined by the follow equation:
%   epsilon = min(epsilon0, K*Vmax)
%               where Vmax is the maximum flux in FVA of the queried
%               direction
%
% USAGE:
%   [epsilon_f, epsilon_r] = makeEpsilonSeq(model, reactions, epsilon0, k)
%
% INPUTS:
%   model:          COBRA model struct
%   reactions:      the list of reactions to calculate epsilons for
%   epsilon0:       the default epsilon (baseline espilon)
%   k:              the coefficient for indicating the custom espilon calling
%                   threshold
%
% OUTPUTS:
%   epsilon_f:      the epsilon sequence for the forward direction of all
%                   input reactions
%   epsilon_r:      the epsilon sequence for the reverse direction of all
%                   input reactions
%   capacity_f:     A logic sequence to indicate whether the forward
%                   direction of the input reaction could carry flux
%   capacity_r:     A logic sequence to indicate whether the reverse
%                   direction of the input reaction could carry flux
% ..AUTHOR: Xuhang Li, 2018

epsilon_f = epsilon0 * ones(length(reactions),1);
epsilon_r = epsilon0 * ones(length(reactions),1);
capacity_f = zeros(length(reactions),1);
capacity_r = zeros(length(reactions),1);

tol = 1e-7; % the default flux numerical tolerance for capacity output
for i = 1: length(reactions)
    testModel = model;
    testModel = changeObjective(testModel,reactions(i));
    Vm_f = optimizeCbModel(testModel, 'max');
    if Vm_f.obj*k < epsilon0 && Vm_f.obj*k > 1e-10 %avoid numeric error
        epsilon_f(i) = Vm_f.obj * k;
    end
    if Vm_f.obj > tol
        capacity_f(i) = true;
    else
        capacity_f(i) = false;
    end

    Vm_r = optimizeCbModel(testModel,'min');
    if -Vm_r.obj*k < epsilon0 && Vm_r.obj*k < -1e-10
        epsilon_r(i) = -Vm_r.obj * k;
    end
    if Vm_r.obj < -tol
        capacity_r(i) = true;
    else
        capacity_r(i) = false;
    end
end
