function solution = optimizeCbModel_dry(model, osenseStr, minNorm, allowLoops, zeroNormApprox)
% this is a dry-run. It returns ready-to-go model structure for gurobi but
% doesn't send it to gurobi API. Will be useful for batch solving!

% Solves a flux balance analysis problem
%
% Solves LP problems of the form
%
% .. math::
%
%    max/min  ~& c^T v \\
%    s.t.     ~& S v = dxdt ~~~~~~~~~~~:y \\
%             ~& C v \leq d~~~~~~~~:y \\
%             ~& lb \leq v \leq ub~~~~:w
%
% USAGE:
%
%    solution = optimizeCbModel(model, osenseStr, minNorm, allowLoops, zeroNormApprox)
%
% INPUT:
%    model:             (the following fields are required - others can be supplied)
%
%                         * S  - `m x 1` Stoichiometric matrix
%                         * c  - `n x 1` Linear objective coefficients
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%
% OPTIONAL INPUTS:
%    model:             
%                         * dxdt - `m x 1` change in concentration with time
%                         * csense - `m x 1` character array with entries in {L,E,G} 
%                           (The code is backward compatible with an m + k x 1 csense vector,
%                           where k is the number of coupling constraints)
%
%                         * C - `k x n` Left hand side of C*v <= d
%                         * d - `k x n` Right hand side of C*v <= d
%                         * dsense - `k x 1` character array with entries in {L,E,G}
%
%    osenseStr:         Maximize ('max')/minimize ('min') (opt, default = 'max')
%    minNorm:           {(0), 'one', 'zero', > 0 , n x 1 vector}, where `[m,n]=size(S)`;
%                       0 - Default, normal LP
%                       'one'  Minimise the Taxicab Norm using LP.
%
%                       .. math::
%
%                          min  ~& |v| \\
%                          s.t. ~& S v = dxdt \\
%                               ~& c^T v = f \\
%                               ~& lb \leq v \leq ub
%
%                       A LP solver is required.
%                       'zero' Minimize the cardinality (zero-norm) of v
%
%                       .. math::
%
%                          min  ~& ||v||_0 \\
%                          s.t. ~& S v = dxdt \\
%                               ~& c^T v = f \\
%                               ~& lb \leq v \leq ub
%
%                       The zero-norm is approximated by a non-convex approximation
%                       Six approximations are available: capped-L1 norm, exponential function
%                       logarithmic function, SCAD function, L_p norm with p<0, L_p norm with 0<p<1
%                       Note : capped-L1, exponential and logarithmic function often give
%                       the best result in term of sparsity.
%
%                       .. See "Le Thi et al., DC approximation approaches for sparse optimization,
%                          European Journal of Operational Research, 2014"
%                          http://dx.doi.org/10.1016/j.ejor.2014.11.031
%                          A LP solver is required.
%
%                       The remaining options work only with a valid QP solver:
%
%                       > 0    Minimises the squared Euclidean Norm of internal fluxes.
%                       Typically 1e-6 works well.
%
%                       .. math::
%
%                          min  ~& 1/2 v'*v \\
%                          s.t. ~& S v = dxdt \\
%                               ~& c^T v = f \\
%                               ~& lb \leq v \leq ub
%
%                       `n` x 1   Forms the diagonal of positive definiate
%                       matrix `F` in the quadratic program
%
%                       .. math::
%
%                          min  ~& 0.5 v^T F v \\
%                          s.t. ~& S v = dxdt \\
%                               ~& c^T v = f \\
%                               ~& lb \leq v \leq ub
%
%    allowLoops:        {0,(1)} If false, then instead of a conventional FBA,
%                       the solver will run an MILP version which does not allow
%                       loops in the final solution.  Default is true.
%                       Runs much slower when set to false.
%                       See `addLoopLawConstraints.m` to for more info.
%
%    zeroNormApprox:    appoximation type of zero-norm (only available when minNorm='zero') (default = 'cappedL1')
%
%                          * 'cappedL1' : Capped-L1 norm
%                          * 'exp'      : Exponential function
%                          * 'log'      : Logarithmic function
%                          * 'SCAD'     : SCAD function
%                          * 'lp-'      : L_p norm with p<0
%                          * 'lp+'      : L_p norm with 0<p<1
%                          * 'all'      : try all approximations and return the best result
%
% OUTPUT:
%    solution:       solution object:
%
%                          * f - Objective value
%                          * v - Reaction rates (Optimal primal variable, legacy FBAsolution.x)
%                          * y - Dual for the metabolites
%                          * w - Reduced costs of the reactions
%                          * s - Slacks of the metabolites
%                          * stat - Solver status in standardized form:
%
%                            * `-1` - No solution reported (timelimit, numerical problem etc)
%                            * `1` - Optimal solution
%                            * `2` - Unbounded solution
%                            * `0` - Infeasible
%                          * origStat - Original status returned by the specific solver
%                    
%                    If the input model contains `C` the following fields are added to the solution:
%
%                          * ctrs_y - the duals for the constraints from C
%                          * ctrs_slack - Slacks of the additional constraints
%
%                    If the model contains the `E` field, the following fields are added to the solution:
%
%                          * vars_v - The optimal primal values of the variables
%                          * vars_w - The reduced costs of the additional variables from E 
%
% .. Author:
%       - Markus Herrgard       9/16/03
%       - Ronan Fleming         4/25/09  Option to minimises the Euclidean Norm of internal
%                                        fluxes using 'cplex_direct' solver
%       - Ronan Fleming         7/27/09  Return an error if any imputs are NaN
%       - Ronan Fleming         10/24/09 Fixed 'E' for all equality constraints
%       - Jan Schellenberger             MILP option to remove flux around loops
%       - Ronan Fleming         12/07/09 Reworked minNorm parameter option to allow
%                                        the full range of approaches for getting
%                                        rid of net flux around loops.
%       - Jan Schellenberger    2/3/09   fixed bug with .f being set incorrectly
%                                        when minNorm was set.
%       - Nathan Lewis          12/2/10  Modified code to allow for inequality
%                                        constraints.
%       - Ronan Fleming         12/03/10 Minor changes to the internal handling of
%                                        global parameters.
%       - Ronan Fleming         14/09/11 Fixed bug in minNorm with negative
%                                        coefficient in objective
%       - Minh Le               11/02/16 Option to minimise the cardinality of
%                                        fluxes vector
%       - Stefania Magnusdottir 06/02/17 Replace LPproblem2 upper bound 10000 with Inf
%       - Ronan Fleming         13/06/17 Support for coupling C*v<=d
%
% NOTE:
%
%    `solution.stat` is either 1, 2, 0 or -1, and is a translation from `solution.origStat`,
%    which is returned by each solver in a solver specific way. That is, not all solvers return
%    the same type of `solution.origStat` and because the cobra toolbox can use many solvers,
%    we need to return to the user of `optimizeCbModel.m` a standard representation, which is what
%    `solution.stat` is.
%
%    When running `optimizeCbModel.m`, unless `solution.stat = 1`, then no solution is returned.
%    This means that it is up to the person calling `optimizeCbModel` to adapt their code to the
%    case when no solution is returned, by checking the value of `solution.stat` first.

if exist('osenseStr', 'var') % Process arguments and set up problem
    if isempty(osenseStr)
        model.osenseStr = 'max';
    else
        model.osenseStr = osenseStr;
    end
else
    if isfield(model, 'osenseStr')
        model.osenseStr = model.osenseStr;
    else
        model.osenseStr = 'max';
    end
end
% Figure out objective sense

if exist('minNorm', 'var')
    if isempty(minNorm)
        %use global solver parameter for minNorm
        minNorm = getCobraSolverParams('LP','minNorm');
    end
    % if minNorm = 'zero' then check the parameter 'zeroNormApprox'
    if isequal(minNorm,'zero')
        if exist('zeroNormApprox', 'var')
            availableApprox = {'cappedL1','exp','log','SCAD','lp-','lp+','all'};
            if ~ismember(zeroNormApprox,availableApprox)
                warning('Approximation is not valid. Default value will be used');
                zeroNormApprox = 'cappedL1';
            end
        else
            zeroNormApprox = 'cappedL1';
        end
    end
else
    %use global solver parameter for minNorm
    minNorm = getCobraSolverParams('LP','minNorm');
end
if exist('allowLoops', 'var')
    if isempty(allowLoops)
        allowLoops = true;
    end
else
    allowLoops = true;
end

%use global solver parameter for printLevel
[printLevel,primalOnlyFlag] = getCobraSolverParams('LP',{'printLevel','primalOnly'});

% size of the stoichiometric matrix
[nMets,nRxns] = size(model.S);

LPproblem = buildLPproblemFromModel(model);

%Double check that all inputs are valid:
if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
    warning('invalid problem');
    return;
end

if isfield(model,'C')
    nCtrs = size(model.C,1);
end

if isfield(model,'E')
    nVars = size(model.E,2);
end

%%
t1 = clock;
% Solve initial LP
if allowLoops
    solution = solveCobraLP_dry(LPproblem);
else
    MILPproblem = addLoopLawConstraints(LPproblem, model, 1:nRxns);
    solution = solveCobraMILP(MILPproblem);
end






