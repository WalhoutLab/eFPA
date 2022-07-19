function [level_f,level_r] = FAA_MILP(MILPproblem_minFlux,model, OFD, epsilon_f,epsilon_r,targetRxns, parforFlag)
% perform FVA given a setup MILProblem and target reactions. The first nRxn
% variables in the MILProblem must be the same as reactions in the input
% model.
%
% USAGE:
%
%    [lb, ub] = FVA_MILP(MILProblem, model, targetRxns,parforFlag)
%
% INPUTS:
%    MILPproblem_minFlux: the input MILP problem (COBRA MILP structure).The
%                       MILP should be readily constrained for FVA
%                       calculation. For example, the total flux cap should
%                       be already set.
%    model:             input model (COBRA model structure)
%    targetRxns:        cell of target reactions to perform FVA on
%    parforFlag:        (0 or 1) whether to use parallel computing
%    BigModel:          (0 or 1) to indicate if the "big model" mode is
%                       used. This mode is recommanded for all complex models. 
%                       In this mode, we release the MILP strigency 
%                       to gain computational speed. But in general, this mode gives almost
%                       identical flux prediction as the normal mode.

%
% OUTPUT:
%   ub:                 a vector of upper boundaries of queried reactions
%   lb:                 a vector of lower boundaries of queried reactions
%
% Additional Notice:    Please make sure the S matrix of the input MILP follows the structure of iMAT++ MILP. Some variables such as absolute flux proxy will be assumed to be at specifc positions, so errors will occur if the S matrix is not formed as standard iMAT++. 
%
% `Yilmaz et al. (2020). Final Tittle and journal.
% .. Author: - Xuhang Li, Mar 2020
if nargin < 6 || isempty(targetRxns)
    targetRxns = model.rxns;
end
if nargin < 7 || isempty(parforFlag)
    parforFlag = true;
end
% if (nargin < 8) 
%     BigModel = 0; % by default, do normal FVA
% end

% if BigModel % This essentially doesnt matter for FAA 
%     relMipGapTol = 1e-4; % we release the MIPgap to 0.1% 
%     % this released MipGap only applies to latent step. The strigency of PFD is still kept.
%     % users can release the MipGap for PFD manually if needed
% else
%     relMipGapTol = 1e-4;
% end
relMipGapTol = 1e-12;
fprintf('Start to perform the FVA...\n');
% Check if is running on gurobi solver
solverOK = changeCobraSolver('gurobi', 'MILP',0);
if ~solverOK
    fprintf('The solver parameter auto-tuning is not supported for current solver! Please use Gurobi for best performance!\n')
end
%% analyze FVA
if parforFlag
    environment = getEnvironment();
    parfor i = 1:length(targetRxns)
        restoreEnvironment(environment);
        targetRxn = targetRxns{i};
        targetInd = strcmp(model.rxns,targetRxn);
        [val_f,type_f, val_r, type_r] = getMinVal(OFD(targetInd), epsilon_f(targetInd),epsilon_r(targetInd));
        % forward direction
        if strcmp(type_f,'NoNeed')
            level_f(i) = 1;
            fprintf('forward direction of %s is in OFD, no need to calculate! \n',targetRxn);
        elseif model.ub(targetInd) <=0
            level_f(i) = -1; % irreversible reaction
            fprintf('forward direction of %s is exluded by boudary, no need to calculate! \n',targetRxn);
        else
            level_f(i) = solveFAA(MILPproblem_minFlux,model,targetRxn, val_f,type_f,'f',solverOK,relMipGapTol);
            fprintf('forward direction of %s found to be of level %d. \n',targetRxn,level_f(i));
        end
        
        % reverse direction
        if strcmp(type_r,'NoNeed')
            level_r(i) = 1;
            fprintf('reverse direction of %s is in OFD, no need to calculate! \n',targetRxn);
        elseif model.lb(targetInd) >=0
            level_r(i) = -1; % irreversible reaction
            fprintf('reverse direction of %s is exluded by boudary, no need to calculate! \n',targetRxn);
        else
            level_r(i) = solveFAA(MILPproblem_minFlux,model,targetRxn, val_r,type_r,'r',solverOK,relMipGapTol);
            fprintf('reverse direction of %s found to be of level %d. \n',targetRxn,level_r(i));
        end
    end
else %same thing but in for loop
    for i = 1:length(targetRxns)
        targetRxn = targetRxns{i};
        targetInd = strcmp(model.rxns,targetRxn);
        [val_f,type_f, val_r, type_r] = getMinVal(OFD(targetInd), epsilon_f(targetInd),epsilon_r(targetInd));
        % forward direction
        if strcmp(type_f,'NoNeed')
            level_f(i) = 1;
         elseif model.ub(targetInd) <=0
            level_f(i) = -1; % irreversible reaction
            fprintf('forward direction of %s is exluded by boudary, no need to calculate! \n',targetRxn);
        else
            level_f(i) = solveFAA(MILPproblem_minFlux,model,targetRxn, val_f,type_f,'f',solverOK,relMipGapTol);
            fprintf('forward direction of %s found to be of level %d. \n',targetRxn,level_f(i));
        end
        
        % reverse direction
        if strcmp(type_r,'NoNeed')
            level_r(i) = 1;
        elseif model.lb(targetInd) >=0
            level_r(i) = -1; % irreversible reaction
            fprintf('reverse direction of %s is exluded by boudary, no need to calculate! \n',targetRxn);
        else
            level_r(i) = solveFAA(MILPproblem_minFlux,model,targetRxn, val_r,type_r,'r',solverOK,relMipGapTol);
            fprintf('reverse direction of %s found to be of level %d. \n',targetRxn,level_r(i));
        end
    end
end
end