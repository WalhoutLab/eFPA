function [status] = solveFAA(MILP,model,targetRxn, val,type,fw,solverOK,relMipGapTol)
fluxTol = 1e-7;
% prepare the expanded A matrix
A = [MILP.A sparse(size(MILP.A,1),1);...
    sparse(1,size(MILP.A,2)) 0];

lbVal = MILP.lb(strcmp(model.rxns,targetRxn));
ubVal = MILP.ub(strcmp(model.rxns,targetRxn));
status = NaN;

if strcmp(fw,'f') % forward direction
    if strcmp(type,'smallFlux') % put rxns that have small fluxes and have to carry flux to OFD
        MILP_tmp = MILP;
        A_tmp = A;
        A_tmp(end,strcmp(model.rxns,targetRxn)) = 1;
        A_tmp(end,end) = fluxTol - ubVal;
        MILP_tmp.A = A_tmp;
        MILP_tmp.csense = [MILP_tmp.csense 'L'];
        MILP_tmp.lb = [MILP_tmp.lb;0];
        MILP_tmp.ub = [MILP_tmp.ub;1];
        MILP_tmp.b = [MILP_tmp.b;fluxTol];
        MILP_tmp.vartype = [MILP_tmp.vartype;'B'];
        MILP_tmp.osense = 1;
        MILP_tmp.c = zeros(size(A_tmp,2),1);
        MILP_tmp.c(end) = 1;
        if ~isempty(MILP_tmp.x0)
            MILP_tmp.x0(end+1) = 1;
        end
        [solution] = autoTuneSolveMILP(MILP_tmp,solverOK,relMipGapTol,targetRxn);
        if ~isnan(solution.obj)
           if  solution.obj >= 1 % I cannot be 0 in the minimization
               status = 1; % qualified, change the status to OFD
           end
        else
            % infeasible error
            % so we are uncertain about the status, go to next step (check
            % for ALT)
        end
    end
    
    if isnan(status) % status not assigned yet (so not a smallflux case or small flux failed to assign)
        % check if in ALT
        MILP_tmp = MILP;
        A_tmp = A;
        A_tmp(end,strcmp(model.rxns,targetRxn)) = 1;
        A_tmp(end,end) = lbVal-val;
        MILP_tmp.A = A_tmp;
        MILP_tmp.csense = [MILP_tmp.csense 'G'];
        MILP_tmp.lb = [MILP_tmp.lb;0];
        MILP_tmp.ub = [MILP_tmp.ub;1];
        MILP_tmp.b = [MILP_tmp.b;lbVal];
        MILP_tmp.vartype = [MILP_tmp.vartype;'B'];
        MILP_tmp.osense = -1;
        MILP_tmp.c = zeros(size(A_tmp,2),1);
        MILP_tmp.c(end) = 1;
        if ~isempty(MILP_tmp.x0)
            MILP_tmp.x0(end+1) = 0;
        end
        [solution] = autoTuneSolveMILP(MILP_tmp,solverOK,relMipGapTol,targetRxn);
        if ~isnan(solution.obj)
           if  solution.obj >= 1
               status = 0;
           else
               status = -1;
           end
        else
            status = 0; % don't block infeasible
        end
    end
    
else % reverse direction
    if strcmp(type,'smallFlux') % rescue rxns that have small fluxes and have to carry flux      
        MILP_tmp = MILP;
        A_tmp = A;
        A_tmp(end,strcmp(model.rxns,targetRxn)) = 1;
        A_tmp(end,end) = -fluxTol - lbVal;
        MILP_tmp.A = A_tmp;
        MILP_tmp.csense = [MILP_tmp.csense 'G'];
        MILP_tmp.lb = [MILP_tmp.lb;0];
        MILP_tmp.ub = [MILP_tmp.ub;1];
        MILP_tmp.b = [MILP_tmp.b;-fluxTol];
        MILP_tmp.vartype = [MILP_tmp.vartype;'B'];
        MILP_tmp.osense = 1;
        MILP_tmp.c = zeros(size(A_tmp,2),1);
        MILP_tmp.c(end) = 1;
        if ~isempty(MILP_tmp.x0)
            MILP_tmp.x0(end+1) = 1;
        end
        [solution] = autoTuneSolveMILP(MILP_tmp,solverOK,relMipGapTol,targetRxn);
        if ~isnan(solution.obj)
           if  solution.obj >= 1 % I cannot be 0 in the minimization
               status = 1; % qualified, change the status
           end
        else
            % infeasible error
            % so we are uncertain about the status, go to next step (check
            % for ALT)
        end
    end
    
    if isnan(status) % status not assigned yet (so not a smallflux case or small flux failed to assign)
        MILP_tmp = MILP;
        A_tmp = A;
        A_tmp(end,strcmp(model.rxns,targetRxn)) = 1;
        A_tmp(end,end) = ubVal+val;
        MILP_tmp.A = A_tmp;
        MILP_tmp.csense = [MILP_tmp.csense 'L'];
        MILP_tmp.lb = [MILP_tmp.lb;0];
        MILP_tmp.ub = [MILP_tmp.ub;1];
        MILP_tmp.b = [MILP_tmp.b;ubVal];
        MILP_tmp.vartype = [MILP_tmp.vartype;'B'];
        MILP_tmp.osense = -1;
        MILP_tmp.c = zeros(size(A_tmp,2),1);
        MILP_tmp.c(end) = 1;
        if ~isempty(MILP_tmp.x0)
            MILP_tmp.x0(end+1) = 0;
        end
        [solution] = autoTuneSolveMILP(MILP_tmp,solverOK,relMipGapTol,targetRxn);
        if ~isnan(solution.obj)
           if  solution.obj >= 1
               status = 0;
           else
               status = -1;
           end
        else
            status = 0; % don't block infeasible
        end
    end  
end
end