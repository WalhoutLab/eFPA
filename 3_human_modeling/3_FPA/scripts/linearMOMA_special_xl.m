function [solutionDel, solStatus] = linearMOMA_special_xl(modelDel, V_wt,MOMAtarget)
% special linear MOMA algrithm for input model who has internal constraints
% also could define the MOMA target
% in the S matrix

[nMets2,~] = size(modelDel.S);
nRxns2 = length(MOMAtarget);
solutionDel.f = [];
solutionDel.x = [];
solutionDel.stat = -1;

% Variables in the following problem are
% x = [v1;v2;delta+;delta-]
% where v1 = wild type flux vector
%       v2 = deletion strain flux vector
%       delta+ = v1 - v2
%       delta- = v2 - v1


% Construct the LHS matrix
% Rows:
% 1: Swt*v1 = 0 for the wild type
% 2: Sdel*v2 = 0 for the deletion strain
% 3: delta+ >= v1-v2
% 4: delta- >= v2-v1
% 5: c'v1 = f1 (wild type)
targetMat = zeros(length(MOMAtarget),length(modelDel.rxns));
for i = 1:length(MOMAtarget)
    targetMat(i,strcmp(MOMAtarget{i},modelDel.rxns)) = 1;
end
A =  [modelDel.S sparse(nMets2,2*nRxns2);
     targetMat eye(nRxns2) sparse(nRxns2,nRxns2);
     -targetMat sparse(nRxns2,nRxns2) eye(nRxns2)];
% Construct the RHS vector
b = [modelDel.b;V_wt; -V_wt];
% Construct the objective (sum of all delta+ and delta-)
c = [zeros(length(modelDel.rxns),1);ones(2*nRxns2,1)];
% Construct the ub/lb
% delta+ and delta- are in [0 100]
lb = [modelDel.lb;zeros(2*nRxns2,1)];
ub = [modelDel.ub;1000*ones(2*nRxns2,1)];
% Construct the constraint direction vector (G for delta's, E for
% everything else)
csense(1:(nMets2)) = modelDel.csense;
csense((nMets2)+1:(nMets2+2*nRxns2)) = 'G';
% Solve the linearMOMA problem
[LPproblem.A,LPproblem.b,LPproblem.c,LPproblem.lb,LPproblem.ub,LPproblem.csense,LPproblem.osense] = deal(A,b,c,lb,ub,csense,1);
LPsolution = solveCobraLP(LPproblem);

if (LPsolution.stat > 0)
    solutionDel.x = LPsolution.full(1:length(modelDel.rxns));
    solutionDel.stat = LPsolution.stat;
    solutionDel.f = LPsolution.obj;
    solutionDel.full = LPsolution.full;
end
solStatus = LPsolution.stat;
