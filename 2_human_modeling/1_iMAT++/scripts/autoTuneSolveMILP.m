function [solution] = autoTuneSolveMILP(MILP,solverOK,relMipGapTol,marker, limitCore)
% The automatic solver optimizer for parallel MILP computation in FVA.
% This function is only designed to minimize the 'false' infeasible
% solution from Gurobi optimizer, and obtain best performance on Gurobi optimizer 
% when doing MILP-based Flux Variability Analysis (FVA) by parfor function in MatLab.
% Other solver or situation may not be applied.
%
% USAGE:
%
%    [solution] = autoTuneSolveMILP(MILP,solverOK,relMipGapTol,marker)
%
% INPUTS:
%    MILP:             	The MILP to be solved (in COBRA format)
%    solverOK:          (0 or 1) indicating if GUROBI solver is avaible
%    relMipGapTol:      the relative MIP gap tolerance
%    marker:            when error occurs, what marker (a string) is
%                       desired to included in the error information print
%
% OUTPUT:
%   solution:           the solution to the input MILP
%
% `Yilmaz, L. S., Li, X., Nanda, S., Fox, B., Schroeder, F., & Walhout, A. J. (2020). Modeling tissue‐relevant Caenorhabditis elegans metabolism at network, pathway, reaction, and metabolite levels. Molecular Systems Biology, 16(10).
% .. Author: - Xuhang Li, June 2020


% non-optimal solution returns NaN in the obj field.

if (nargin < 5) 
    limitCore = 1; % by default, we assume parfor is used and we limit the solver core usage
end

if solverOK && limitCore% for best parfor parallel performance, we limit 1 core per GUROBI thread 
    gurobiParameters = struct();
    gurobiParameters.Threads = 1;
    solution = solveCobraMILP_XL(MILP, gurobiParameters,'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
else % go with default for other solvers
    solution = solveCobraMILP_XL(MILP, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
end

if solution.stat ==0 && solverOK% when failed to solve (infeasible), we start to tune solver parameter 
    %NOTE: SPECIFIC TO GUROBI SOLVER!
    gurobiParameters = struct();
    gurobiParameters.Presolve = 0;
    % we dont limit to 1 core to push solver to search solution space
    solution = solveCobraMILP_XL(MILP,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
    if solution.stat == 0 % still infeasible, we force maximum numeric focus
        gurobiParameters.NumericFocus = 3;
        solution = solveCobraMILP_XL(MILP,gurobiParameters, 'timeLimit', 1200, 'logFile', 'MILPlog', 'printLevel', 0,'relMipGapTol',relMipGapTol);
    end
end

if solution.stat ~= 1 % give up and report error
    fprintf('infeasible or unbounded occured!(%s) (error code: %d)\n',marker,solution.stat);
    solution.obj = nan;
end
end
