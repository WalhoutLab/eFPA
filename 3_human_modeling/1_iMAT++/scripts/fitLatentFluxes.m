function [FluxDistribution, latentRxn,Nfit_latent, minTotal,MILPproblem_minFlux] = fitLatentFluxes(MILProblem, model, PFD, Hgenes, epsilon_f,epsilon_r,latentCAP,verbose,relMipGapTol,sideCAP)
% The latent flux fitting module in iMAT++ pipeline. This function
% works with a specific COBRA MILP input that is the product of primary
% iMAT++ fitting. This function will find latent reactions and apply
% fitting procedure similar to iMAT algorithm while only maximizing the
% number of latent reactions that carry flux. This function is not a
% standalone integration function. Please see "IMATplusplus" for details.
%
% USAGE:
%
%    FluxDistribution = fitLatentFluxes(MILProblem, model, PFD, Hgenes, epsilon_f,epsilon_r,latentCAP)
%
% INPUTS:
%    MILProblem:        the input MILP problem (COBRA MILP structure). The
%                       input MILP should be already constrained with Nfit,
%                       MinLow and total flux cap; Besides, the last row of
%                       the S matrix must be the constraint for total flux
%                       cap, since we will iteratively modify it in the
%                       latent fitting procedure.
%    model:             input model (COBRA model structure)
%    PFD:               Primary Flux Distribution 
%    Hgenes:            the list of highly expressed genes
%    epsilon_f:         the epsilon sequence for the forward direction of all reactions. Non-applicable reactions (i.e., for a irreversible,
%                       reverse-only reaction) should have a non-zero default to avoid numeric error.
%    epsilon_r:         the epsilon sequence for the reverse direction of all reactions. Non-applicable reactions (i.e., for a irreversible,
%                       forward-only reaction) should have a non-zero default to avoid numeric error.
%    latentCAP:         the total flux cap for recursive fitting of latent
%                       reactions. The total flux will be capped at (1 +
%                       latentCAP)*OriginalTotalFlux; The default cap is
%                       0.01
%    verbose:           (0 or 1) show the MILP solving details or not
%
% OUTPUT:
%   FluxDistribution:   the resulting flux distribution, Optimal Flux Distribution (OFD)
% OPTIONAL OUTPUTS:
%   latentRxn:          the list of all identified latent reactions to be fitted 
%   Nfit_latent:        the (total) number of latent reactions fitted
%   minTotal:           the minimal total flux of OFD (Optimal Flux Distribution)
%   MILPproblem_minFlux:the final MILP structure in the last latent calculation
%
% Additional Notice:    Please make sure the S matrix of the input MILP follows the structure of iMAT++ MILP format. Some variables such as absolute flux proxy will be assumed to be at specifc positions, so errors will occur if the S matrix is not formed as standard iMAT++. 
%
% `Yilmaz, L. S., Li, X., Nanda, S., Fox, B., Schroeder, F., & Walhout, A. J. (2020). Modeling tissue‐relevant Caenorhabditis elegans metabolism at network, pathway, reaction, and metabolite levels. Molecular Systems Biology, 16(10).
% .. Author: - Xuhang Li, Mar 2020
if nargin < 8 || isempty(verbose)
    verbose = 0;
end
% Check if is running on gurobi solver
solverOK = changeCobraSolver('gurobi', 'MILP',0);
if ~solverOK
    fprintf('The solver parameter auto-tuning is not supported for the current solver! Please use Gurobi for best performance!\n')
end
%% store the constraint index for total flux
minTotalInd = length(MILProblem.b); % assume the total flux constriant is the last row!
fprintf('the total flux is constrained to %.2f \n',MILProblem.b(end));
%% define candidate reactions of latent reactions from HGene list
latentCandi = {};
for i = 1:length(model.rxns)
    mygenes = model.genes(logical(model.rxnGeneMat(i,:)));
    if all(ismember(mygenes,Hgenes)) && ~isempty(mygenes) % note that empty is not desired!
        latentCandi = [latentCandi;model.rxns(i)];
    end
end
fprintf('starting the latent fitting loop... \n');
%% start of the latent fitting loop
latentRxn = {};
count = 0;
while 1
    tic()
    %% define latent reactions
    % get potentially flux carrying reaction list 
    fluxProduct = model.S .* PFD';
    MetFlux = sum(abs(fluxProduct),2) / 2; % total flux for each metabolite
    FluxCapMet = model.mets(MetFlux >= 1e-5); % the numerical tolerance is 1e-5; these met carries valid flux
    % get the new latent rxns set
    newLatent = {};
    for i = 1:length(latentCandi)
        metsF = model.mets(model.S(:,strcmp(model.rxns,latentCandi{i}))<0);
        metsF(cellfun(@(x) ~isempty(regexp(x,'NonMetConst', 'once')),metsF)) = [];% for tissue model only, but won't cause error for other model. We remove the special mets (which are place holder for constraints)
        metsR = model.mets(model.S(:,strcmp(model.rxns,latentCandi{i}))>0);
        metsR(cellfun(@(x) ~isempty(regexp(x,'NonMetConst', 'once')),metsR)) = [];
        % we require either side of a reaction to have active metabolite to
        % call a latent reaction
        if ((model.lb(strcmp(model.rxns,latentCandi{i}))<0 && model.ub(strcmp(model.rxns,latentCandi{i}))>0) && ((all(ismember(metsF,FluxCapMet))&& ~isempty(metsF)) || (all(ismember(metsR,FluxCapMet))&& ~isempty(metsR))))...% reversible & one side mets all have flux
           ||...
            ((model.lb(strcmp(model.rxns,latentCandi{i}))>=0) && (all(ismember(metsF,FluxCapMet))&& ~isempty(metsF)))...% irreversible & reactant side mets all have flux    
           ||...
            ((model.ub(strcmp(model.rxns,latentCandi{i}))<=0) && (all(ismember(metsR,FluxCapMet))&& ~isempty(metsR))) % irreversible & product side mets all have flux
            newLatent = [newLatent;latentCandi(i)];
        end
    end
    newLatent = newLatent(~ismember(newLatent,model.rxns(abs(PFD)>=1e-5)));% only the closed reaction in PFD are latent
    % union new latent reaction with previous latent reaction!
    % no new latent rxns? ==> stop the iteration
    if isempty(setdiff(newLatent,latentRxn))
        fprintf('...done! \n');
        break;
    else
        fprintf('...found %d new latent rxns \n',length(setdiff(newLatent,latentRxn)));
        latentRxn = union(latentRxn,newLatent);
    end
    %% maximaze number of "active" latent reaction (like iMAT but without forcing zero flux to lowly expressed reactions!)
    [A B] = ismember(latentRxn,model.rxns);
    latentInd = B(A);
    epsilon_f_sorted = epsilon_f(latentInd);
    epsilon_r_sorted = epsilon_r(latentInd);
    % Creating A matrix
    A = [MILProblem.A sparse(size(MILProblem.A,1),2*length(latentRxn));...
        sparse(2*length(latentRxn),size(MILProblem.A,2)) sparse(2*length(latentRxn),2*length(latentRxn))];
    for i = 1:length(latentRxn)
        A(i+size(MILProblem.A,1),latentInd(i)) = 1;
        A(i+size(MILProblem.A,1),i+size(MILProblem.A,2)) = model.lb(latentInd(i)) - epsilon_f_sorted(i);
        A(i+size(MILProblem.A,1)+length(latentRxn),latentInd(i)) = 1;
        A(i+size(MILProblem.A,1)+length(latentRxn),i+size(MILProblem.A,2)+length(latentRxn)) = model.ub(latentInd(i)) + epsilon_r_sorted(i);
    end
    % variable type
    vartype1(1:size(MILProblem.A,2),1) = MILProblem.vartype;
    vartype2(1:2*length(latentRxn),1) = 'B';
    vartype = [vartype1;vartype2];
    % Creating csense
    csense1(1:size(MILProblem.A,1)) = MILProblem.csense;
    csense2(1:length(latentRxn)) = 'G';
    csense3(1:length(latentRxn)) = 'L';
    csense = [csense1 csense2 csense3];
    % Creating lb and ub
    lb_y = zeros(2*length(latentRxn),1);
    ub_y = ones(2*length(latentRxn),1);
    lb = [MILProblem.lb;lb_y];
    ub = [MILProblem.ub;ub_y];
    % Creating c
    c_v = zeros(size(MILProblem.A,2),1);
    c_y = ones(2*length(latentRxn),1);
    c = [c_v;c_y];
    % Creating b
    b_s = MILProblem.b;
    lb_rh = model.lb(latentInd);
    ub_rh = model.ub(latentInd);
    b = [b_s;lb_rh;ub_rh];

    MILPproblem_latent.A = A;
    MILPproblem_latent.b = b;
    MILPproblem_latent.c = c;
    MILPproblem_latent.lb = lb;
    MILPproblem_latent.ub = ub;
    MILPproblem_latent.csense = csense;
    MILPproblem_latent.vartype = vartype;
    MILPproblem_latent.osense = -1;
    % define the initial solution
    extraX0_1 = (MILProblem.x0(latentInd) >= epsilon_f_sorted)*1;
    extraX0_2 = (MILProblem.x0(latentInd) <= -epsilon_r_sorted)*1;
    MILPproblem_latent.x0 = [MILProblem.x0(1:length(MILProblem.vartype)); extraX0_1; extraX0_2];
    solution = solveCobraMILP_XL(MILPproblem_latent, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
    if solution.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter
        % NOTE: SPECIFIC TO GUROBI SOLVER!
        fprintf('Initial solving failed! Start the auto-tune...#1\n')
        gurobiParameters = struct();
        gurobiParameters.Presolve = 0;
        solution = solveCobraMILP_XL(MILPproblem_latent,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
        if solution.stat == 0
            fprintf('Initial solving failed! Start the auto-tune...#2\n')
            gurobiParameters.NumericFocus = 3;
            solution = solveCobraMILP_XL(MILPproblem_latent,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
            if solution.stat == 0
                error('MILP solving failed! Please inspect the reason!');
            end
        end
    end   
    if solution.stat ~= 1
        if solution.stat == 3
            % pass with warning
            warning('Latent reactions are not fully fitted! (solver time out! Change timeLimit if needed!)');
        else % unknown error
            error('MILP solving failed! Please inspect the reason!');
        end
    end
    fprintf('...total latent fitted: %d \n',sum(solution.int(end-2*length(latentRxn)+1:end)));
    Nfit_latent = sum(solution.int(end-2*length(latentRxn)+1:end));
    %% minimize total flux to update the PFD
    MILPproblem_minFlux = solution2constraint(MILPproblem_latent,solution);
    MILPproblem_minFlux.x0 = solution.full;
    % erase the total flux constraint (otherwise nutrient uptake may not be
    % sufficiently minimized)
    MILPproblem_minFlux.b(minTotalInd) = 0; % we change it to total flux > 0 (always true)
    MILPproblem_minFlux.csense(minTotalInd) = 'G';
    % minimize total flux
    %% 12052020 add major nutrient uptake constriants - for human integration
    if any(strcmp(model.rxns,'EX_sideNutr')) % run the dynamic constriant for unlimited bacUp
        MILProblemTemp = MILPproblem_minFlux;
        MILProblemTemp.c = zeros(length(MILPproblem_minFlux.c),1);
        MILProblemTemp.c(find(strcmp(model.rxns,'EX_sideNutr'))) = 1;
        MILProblemTemp.osense = 1;
        solution_tmp = solveCobraMILP_XL(MILProblemTemp, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);        
        if solution_tmp.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter 
            fprintf('Initial solving failed! Start the auto-tune...#1\n')
            gurobiParameters = struct();
            gurobiParameters.Presolve = 0;
            solution_tmp = solveCobraMILP_XL(MILProblemTemp,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
            if solution_tmp.stat == 0
                fprintf('Initial solving failed! Start the auto-tune...#2\n')
                gurobiParameters.NumericFocus = 3;
                solution_tmp = solveCobraMILP_XL(MILProblemTemp,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
                if solution_tmp.stat == 0
                    error('MILP solving failed! Please inspect the reason!');
                end
            end
        end   
        fprintf('the minimal side nutrient demand is %f \n',solution_tmp.obj);
        % impose boundary constraint
        MILPproblem_minFlux.ub(find(strcmp(model.rxns,'EX_sideNutr'))) = sideCAP * solution_tmp.obj;
        model.ub(find(strcmp(model.rxns,'EX_sideNutr'))) = sideCAP * solution_tmp.obj;
    end
    
    % NOTE: this section is specific to the MILP structure in iMAT++! 
    % We use the absolute flux proxy variables in the INPUT MILProblem instead of creating new variables
    % create a new objective function
    % Creating c vector (objective function)
    c_minFlux = zeros(size(MILProblem.A,2),1);
    c_minFlux(end-length(model.rxns)+1:end) = ones(length(model.rxns),1);
    c = [c_minFlux;zeros(2*length(latentRxn),1)];
    MILPproblem_minFlux.c = c;
    MILPproblem_minFlux.osense = 1;
    solution = solveCobraMILP_XL(MILPproblem_minFlux, 'timeLimit', 1200, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
    if solution.stat == 0 && solverOK% when failed to solve, we start to tune solver parameter 
        fprintf('Initial solving failed! Start the auto-tune...#1\n')
        gurobiParameters = struct();
        gurobiParameters.Presolve = 0;
        solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
        if solution.stat == 0
            fprintf('Initial solving failed! Start the auto-tune...#2\n')
            gurobiParameters.NumericFocus = 3;
            solution = solveCobraMILP_XL(MILPproblem_minFlux,gurobiParameters, 'timeLimit', 600, 'logFile', 'MILPlog', 'printLevel', verbose,'relMipGapTol',relMipGapTol);
            if solution.stat == 0
                error('MILP solving failed! Please inspect the reason!');
            end
        end
    end   
    if solution.stat ~= 1
        if solution.stat == 3
            % pass with warning
            warning('Latent total flux is not fully minimized! (solver time out! Change timeLimit if needed!)');
        else % unknown error
            error('MILP solving failed! Please inspect the reason!');
        end
    end
    minTotal = solution.obj;
    PFD = solution.full(1:length(model.rxns));
    count = count +1;
    fprintf('latent fitting cycle completed .... #%d \n',count);
    %% update the minTotal constraint
    MILProblem.b(minTotalInd) = minTotal*(1+latentCAP);
    MILProblem.x0 = solution.full;% update the initial solution
    % if applicable, also the major nutrient constraint
    if any(strcmp(model.rxns,'EX_sideNutr'))
        MILProblem.ub(find(strcmp(model.rxns,'EX_sideNutr'))) = 1000;
        model.ub(find(strcmp(model.rxns,'EX_sideNutr'))) = 1000;
    end
    
    fprintf('the total flux constraint is updated to %.2f \n',MILProblem.b(minTotalInd));
    toc()
end
% return the OFD (last updated PFD is OFD)
FluxDistribution = PFD;
% return the MILP (the total flux is constrained by the 5% cap of the total flux of OFD)
MILPproblem_minFlux.b(minTotalInd) = minTotal*(1+latentCAP);
MILPproblem_minFlux.x0 = solution.full;% start from minimal flux state may speed up MILP solving
end
