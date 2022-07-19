function solutions = solveGurobi_batch(LPproblem_batch,osense)
% build the disconnected model
LPproblem = struct();
LPproblem.A = sparse(sum(cellfun(@(x) size(x.A,1),LPproblem_batch)),sum(cellfun(@(x) size(x.A,2),LPproblem_batch)));
LPproblem.ub = zeros(sum(cellfun(@(x) length(x.ub),LPproblem_batch)),1);
LPproblem.lb = zeros(sum(cellfun(@(x) length(x.lb),LPproblem_batch)),1);
LPproblem.c = zeros(sum(cellfun(@(x) length(x.c),LPproblem_batch)),1);
LPproblem.b = zeros(sum(cellfun(@(x) length(x.b),LPproblem_batch)),1);
LPproblem.osense = osense;
LPproblem.modelsense = osense;
LPproblem.csense = blanks(sum(cellfun(@(x) length(x.csense),LPproblem_batch)))';
LPproblem.rhs = zeros(sum(cellfun(@(x) length(x.rhs),LPproblem_batch)),1);
LPproblem.obj = zeros(sum(cellfun(@(x) length(x.obj),LPproblem_batch)),1);
LPproblem.sense = blanks(sum(cellfun(@(x) length(x.sense),LPproblem_batch)))';

LPproblem.A(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1))), 1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1))))...
     = LPproblem_batch{1}.A;
LPproblem.ub(1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1)))) = LPproblem_batch{1}.ub;
LPproblem.lb(1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1)))) = LPproblem_batch{1}.lb;
LPproblem.c(1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1)))) = LPproblem_batch{1}.c;
LPproblem.b(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1)))) = LPproblem_batch{1}.b;
LPproblem.csense(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1)))) = LPproblem_batch{1}.csense;
LPproblem.rhs(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1)))) = LPproblem_batch{1}.rhs;
LPproblem.obj(1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1)))) = LPproblem_batch{1}.obj;
LPproblem.sense(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1)))) = LPproblem_batch{1}.sense;
for i = 2:length(LPproblem_batch)
    startx = sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1:i-1)));
    starty = sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1:i-1)));
    lenx = size(LPproblem_batch{i}.A,1);
    leny = size(LPproblem_batch{i}.A,2);
    LPproblem.A(startx+1:startx+lenx, starty+1:starty+leny) = LPproblem_batch{i}.A;
    LPproblem.ub(starty+1:starty+leny) = LPproblem_batch{i}.ub;
    LPproblem.lb(starty+1:starty+leny) = LPproblem_batch{i}.lb;
    LPproblem.c(starty+1:starty+leny) = LPproblem_batch{i}.c;
    LPproblem.b(startx+1:startx+lenx) = LPproblem_batch{i}.b;
    LPproblem.csense(startx+1:startx+lenx) = LPproblem_batch{i}.csense;
    LPproblem.rhs(startx+1:startx+lenx) = LPproblem_batch{i}.rhs;
    LPproblem.obj(starty+1:starty+leny) = LPproblem_batch{i}.obj;
    LPproblem.sense(startx+1:startx+lenx) = LPproblem_batch{i}.sense;
end

% load default parameter
% get the solver parameters
[cobraParams, solverParams] = parseSolverParameters('LP');
param=solverParams;
if ~isfield(param,'FeasibilityTol')
        param.FeasibilityTol = cobraParams.feasTol;
end
if ~isfield(param,'OptimalityTol')
    param.OptimalityTol = cobraParams.optTol;
end
param.OutputFlag = 0;
% set the solver specific parameters
param = updateStructData(param,solverParams);
param.Threads = 1;
param.TimeLimit = length(LPproblem_batch) * 10;
% call the solver
resultgurobi = gurobi(LPproblem,param);

% switch back to numeric
if strcmp(LPproblem.osense,'max')
    LPproblem.osense = -1;
else
    LPproblem.osense = 1;
end
% see the solvers original status -Ronan
origStat = resultgurobi.status;
if strcmp(resultgurobi.status, 'OPTIMAL')
    stat = 1; % optimal solution found
    [x,f,y,w] = deal(resultgurobi.x,resultgurobi.objval,LPproblem.osense*resultgurobi.pi,LPproblem.osense*resultgurobi.rc);

    s = LPproblem.b - LPproblem.A * x; % output the slack variables

    % save the basis
    basis.vbasis=resultgurobi.vbasis;
    basis.cbasis=resultgurobi.cbasis;

    % assign solution
    [solution.full, solution.obj, solution.rcost, solution.dual, solution.slack, ...
     solution.solver, solution.stat, solution.origStat, ...
     solution.basis] = deal(x,f,w,y,s,'Gurobi',stat,origStat,basis);

    %% split solutions
    solutions = {};

    solutions{1} = struct();
    solutions{1}.full = solution.full(1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1))));
    solutions{1}.rcost = solution.rcost(1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1))));
    solutions{1}.dual = solution.dual(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1))));
    solutions{1}.slack= solution.slack(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1))));
    solutions{1}.basis = struct();
    solutions{1}.basis.vbasis = solution.basis.vbasis(1:sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1))));
    solutions{1}.basis.cbasis = solution.basis.vbasis(1:sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1))));
    solutions{1}.solver = 'Gurobi';
    solutions{1}.stat = solution.stat;
    solutions{1}.origStat = solution.origStat;
    solutions{1}.obj = solutions{1}.full' * LPproblem_batch{1}.obj;  
    for i = 2:length(LPproblem_batch)
        starty = sum(cellfun(@(x) size(x.A,1),LPproblem_batch(1:i-1)));
        startx = sum(cellfun(@(x) size(x.A,2),LPproblem_batch(1:i-1)));
        leny = size(LPproblem_batch{i}.A,1);
        lenx = size(LPproblem_batch{i}.A,2);
        solutions{i} = struct();
        solutions{i}.full = solution.full(startx+1:startx+lenx);
        solutions{i}.rcost = solution.rcost(startx+1:startx+lenx);
        solutions{i}.dual = solution.dual(starty+1:starty+leny);
        solutions{i}.slack= solution.slack(starty+1:starty+leny);
        solutions{i}.basis = struct();
        solutions{i}.basis.vbasis = solution.basis.vbasis(startx+1:startx+lenx);
        solutions{i}.basis.cbasis = solution.basis.vbasis(starty+1:starty+leny);
        solutions{i}.solver = 'Gurobi';
        solutions{i}.stat = solution.stat;
        solutions{i}.origStat = solution.origStat;
        solutions{i}.obj = solutions{i}.full' * LPproblem_batch{i}.obj;  
    end
else % lumped problem is numerically hard, so we solve seperately
    solutions = {};
    for i = 1:length(LPproblem_batch)
        LPproblem_batch{i}.osense = osense;
        resultgurobi = gurobi(LPproblem_batch{i},param);
        
        % switch back to numeric
        if strcmp(LPproblem_batch{i}.osense,'max')
            LPproblem_batch{i}.osense = -1;
        else
            LPproblem_batch{i}.osense = 1;
        end
        origStat = resultgurobi.status;
        
        if strcmp(resultgurobi.status, 'OPTIMAL')
            stat = 1; % optimal solution found
            [x,f,y,w] = deal(resultgurobi.x,resultgurobi.objval,LPproblem_batch{i}.osense*resultgurobi.pi,LPproblem_batch{i}.osense*resultgurobi.rc);
            s = LPproblem_batch{i}.b - LPproblem_batch{i}.A * x; % output the slack variables
            % save the basis
            basis.vbasis=resultgurobi.vbasis;
            basis.cbasis=resultgurobi.cbasis;
            % assign solution
            [solutions{i}.full, solutions{i}.obj, solutions{i}.rcost, solutions{i}.dual, solutions{i}.slack, ...
             solutions{i}.solver, solutions{i}.stat, solutions{i}.origStat, ...
             solutions{i}.basis] = deal(x,f,w,y,s,'Gurobi',stat,origStat,basis);
        else % something is wrong
            resultgurobi
            error('FPA failed, check reason!')
        end
    end
end

