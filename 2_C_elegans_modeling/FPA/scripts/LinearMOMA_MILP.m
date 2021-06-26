function [solutionDel, solStat] = LinearMOMA_MILP(modelDel,MILP, V_wt,MOMAtarget)
% special linear MOMA algrithm for a constrained MILP input
% also could define the MOMA target (but note: the v_wt should be aligned
% with MOMA target!!)
% in the S matrix
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

[nMets2,~] = size(modelDel.S);
nRxns2 = length(MOMAtarget);
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
A =  [MILP.A sparse(size(MILP.A,1),2*nRxns2);
     targetMat sparse(nRxns2,size(MILP.A,2)-size(targetMat,2)) eye(nRxns2) sparse(nRxns2,nRxns2);
     -targetMat sparse(nRxns2,size(MILP.A,2)-size(targetMat,2)) sparse(nRxns2,nRxns2) eye(nRxns2)];
% Construct the RHS vector
b = [MILP.b;V_wt; -V_wt];
% Construct the objective (sum of all delta+ and delta-)
c = [zeros(size(MILP.A,2),1);ones(2*nRxns2,1)];
% Construct the ub/lb
% delta+ and delta- are in [0 2000]
lb = [MILP.lb;zeros(2*nRxns2,1)];
ub = [MILP.ub;2000*ones(2*nRxns2,1)];
% Construct the constraint direction vector (G for delta's, E for
% everything else)
csense(1:size(MILP.A,1)) = MILP.csense;
csense((size(MILP.A,1))+1:(size(MILP.A,1)+2*nRxns2)) = 'G';
% construct other input for a MILP
x0 = [];
vartype1(1:2*nRxns2,1) = 'C';
vartype = [MILP.vartype;vartype1];
% Solve the linearMOMA problem
[MILPproblem.A,MILPproblem.b,MILPproblem.c,MILPproblem.lb,MILPproblem.ub,MILPproblem.csense,MILPproblem.osense,MILPproblem.vartype,MILPproblem.x0] = deal(A,b,c,lb,ub,csense,1,vartype,x0);
try 
    MILPsolution = solveCobraMILP_XL(MILPproblem, 'timeLimit', 7200, 'logFile', 'MILPlog', 'printLevel', 0);
    solutionDel.x = MILPsolution.full(1:length(modelDel.rxns));
    solutionDel.stat = MILPsolution.stat;
    solutionDel.obj = MILPsolution.obj;
catch
    solutionDel.stat = -3; %-3 => solver error!
end
solStat = solutionDel.stat;
end
