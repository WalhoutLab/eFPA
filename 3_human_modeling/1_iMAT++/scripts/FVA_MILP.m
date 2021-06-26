function [FVA_lb, FVA_ub] = FVA_MILP(MILPproblem_minFlux, model, targetRxns,parforFlag,relMipGapTol)
% Perform Flux Variability Analysis (FVA) given a MILProblem and target reactions. 
% The first nRxn variables in the MILProblem must be the fluxes of reactions in the input
% model.
%
% USAGE:
%
%    [lb, ub] = FVA_MILP(MILProblem, model, targetRxns, parforFlag)
%
% INPUTS:
%    MILPproblem_minFlux: the input MILP problem (COBRA MILP structure). 
%                       The MILP should be readily constrained for FVA
%                       calculation. For example, the total flux cap should
%                       be already set.
%    model:             input model (COBRA model structure)
%    targetRxns:        cell of target reactions to perform FVA on
%    parforFlag:        (0 or 1) whether to use parallel computing (parfor)
%    relMipGapTol:      the relative MIP gap to use in the FVA calculation
%
% OUTPUT:
%   ub:                 a vector of upper boundaries of queried reactions
%   lb:                 a vector of lower boundaries of queried reactions
%
% Additional Notice:    Please make sure the S matrix of the input MILP follows the structure of iMAT++ MILP. Some variables such as absolute flux proxy will be assumed to be at specifc positions, so errors will occur if the S matrix is not formed as standard iMAT++. 
%
% `Yilmaz, L. S., Li, X., Nanda, S., Fox, B., Schroeder, F., & Walhout, A. J. (2020). Modeling tissue‐relevant Caenorhabditis elegans metabolism at network, pathway, reaction, and metabolite levels. Molecular Systems Biology, 16(10).
% .. Author: - Xuhang Li, Mar 2020
if nargin < 3 || isempty(targetRxns)
    targetRxns = model.rxns;
end
if nargin < 4 || isempty(parforFlag)
    parforFlag = true;
end
if (nargin < 5) 
    relMipGapTol = 1e-12;
end

fprintf('Start to perform the FVA...\n');
% Check if is running on gurobi solver
solverOK = changeCobraSolver('gurobi', 'MILP',0);
if ~solverOK
    fprintf('The solver parameter auto-tuning is not supported for current solver! Please use Gurobi for best performance!\n')
end
%% analyze FVA
if parforFlag
    environment = getEnvironment();
    MILPproblem_minFlux_ori = MILPproblem_minFlux;
    parfor i = 1:length(targetRxns)
        restoreEnvironment(environment);
        MILPproblem_minFlux = MILPproblem_minFlux_ori;
        targetRxn = targetRxns(i);
        FluxObj = find(ismember(model.rxns,targetRxn)); 
        % create a new objective function
        c = zeros(size(MILPproblem_minFlux.A,2),1);
        c(FluxObj) = 1;
        MILPproblem_minFlux.c = c;
        % reverse direction (lb)
        MILPproblem_minFlux.osense = 1;
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn{:});       
        FVA_lb(i) = solution.obj;
        fprintf('lower boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        % forward direction (ub)
        MILPproblem_minFlux.osense = -1;
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn{:}); 
        FVA_ub(i) = solution.obj;
        fprintf('upper boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
    end
else % same thing but in a for loop
    for i = 1:length(targetRxns)
        targetRxn = targetRxns(i);
        FluxObj = find(ismember(model.rxns,targetRxn)); 
        % create a new objective function
        c = zeros(size(MILPproblem_minFlux.A,2),1);
        c(FluxObj) = 1;
        MILPproblem_minFlux.c = c;
        MILPproblem_minFlux.osense = 1;
        % when parfor is not used, go with default (use up cores)
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn{:},0);       
        FVA_lb(i) = solution.obj;
        fprintf('lower boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
        MILPproblem_minFlux.osense = -1;
        % solve the MILP
        solution = autoTuneSolveMILP(MILPproblem_minFlux,solverOK,relMipGapTol,targetRxn{:},0); 
        FVA_ub(i) = solution.obj;
        fprintf('upper boundary of %s found to be %f. \n',targetRxn{:},solution.obj);
    end
end
end
