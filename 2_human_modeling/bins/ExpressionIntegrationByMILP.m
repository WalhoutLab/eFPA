function [OpenGene, OpenedHReaction,ClosedLReaction,solution, MILPproblem] = ExpressionIntegrationByMILP(model, RHNames, RLNames,HGeneList,epsilon_f, epsilon_r, logfile, runtime,verbose)
% the function is based on the original "iMAT.m" function in COBRA toolbox
% with significant alteration. It implements the first step fitting of
% highly and rarely expressed genes in IMAT++ algorithm. In short, it converts the gene expression
% data into a MILP where number of flux-carrying highly expressed genes is
% maximized together with flux-depeleted rarely expressed reactions
%
% USAGE:
%
%    [OpenGene, OpenedHReaction,ClosedLReaction,solution, MILPproblem] = ExpressionIntegrationByMILP(model, RHNames, RLNames,HGeneList,epsilon_f, epsilon_r, logfile, runtime)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    RHNames:           cell with reaction names that are labled as highly expressed
%    RLNames:           cell with reaction names that are labled as rarely
%                       expressed; actually is zero level in the
%                       quadra-level catagory and needs to be minimized
%    HGeneList:         cell with gene names that are labeled as highly
%                       expressed
%    epsilon_f:         the epsilon sequence for the forward direction of
%                       all reactions in the model (needs to be the same length of model.rxns)
%    epsilon_r:         the epsilon sequence for the reverse direction of
%                       all reactions in the model (needs to be the same length of model.rxns)
%
% OPTIONAL INPUTS:
%    logfile:           name of the file to save the MILP log (string)
%    runtime:           maximum solve time for the MILP (default value - 7200s)
%    verbose:           (0 or 1) show the MILP solving details or not
%
% OUTPUT:
%    OpenGene:       	a list of highly expressed genes that are successfully fitted
%                       (carrying flux)
%    OpenedHReaction:   a list of reactions that are both labeled as high
%                       reaction and successfully fitted (carrying flux)
%    ClosedLReaction:    a list of reactions that are both labeled as rarely expressed
%                       reaction and successfully fitted (carrying no flux)
%    solution:          the solution struct of the MILProblem
%    MILPproblem:       the constructed MILProblem for IMAT++ fitting
%
% ..Author:     Xuhang Li, Dec 2018 

if nargin < 7 || isempty(runtime)
    runtime = 300;
end
if nargin < 8 || isempty(logfile)
    logfile = 'MILPlog';
end
if nargin < 9 || isempty(verbose)
    verbose = 0;
end
%% step1: construct the index for highly expressed reactions (RHindex) and rarely expressed rxns
RHindex = find(ismember(model.rxns, RHNames));
RLindex = find(ismember(model.rxns, RLNames));
% only epsilons for Hrxns will be used, and we need to make sure the order
% of epsilons matches Hrxns
epsilon_f_sorted = epsilon_f(ismember(model.rxns,RHNames));
epsilon_r_sorted = epsilon_r(ismember(model.rxns,RHNames));

%% step2: construct the MILP for fitting
% get the inputs for making the MILP
S = model.S;
lb = model.lb;
ub = model.ub;
% construct the basic MILP with objective function undefined
% -------------------------
% Creating A matrix
A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end
% constructing constraints for fitting high reactions
for i = 1:length(RHindex)
    A(i+size(S,1),RHindex(i)) = 1;
    A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - epsilon_f_sorted(i);
    A(i+size(S,1)+length(RHindex),RHindex(i)) = 1;
    A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + epsilon_r_sorted(i);
end
% constructing constraints for fitting zero (rarely expressed) reactions
for i = 1:length(RLindex)
    A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1;
    A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i));
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1;
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i));
end
% Creating csense
csense1(1:size(S,1)) = model.csense;
csense2(1:length(RHindex)) = 'G';
csense3(1:length(RHindex)) = 'L';
csense4(1:length(RLindex)) = 'G';
csense5(1:length(RLindex)) = 'L';
csense = [csense1 csense2 csense3 csense4 csense5];
% Creating lb and ub
lb_y = zeros(2*length(RHindex)+length(RLindex),1);
ub_y = ones(2*length(RHindex)+length(RLindex),1);
lb = [lb;lb_y];
ub = [ub;ub_y];
% Creating b
b_s = model.b;
lb_rh = lb(RHindex);
ub_rh = ub(RHindex);
lb_rl = lb(RLindex);
ub_rl = ub(RLindex);
b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl];
% Creating vartype
vartype1(1:size(S,2),1) = 'C';
vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
vartype = [vartype1;vartype2];
% create MILP
MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.csense = csense;
MILPproblem.vartype = vartype;
MILPproblem.osense = -1;
MILPproblem.x0 = [];
% ---------------------------
%% step3: expand the A matrix to add variables for gene-centric integration
% add the variables (columns in A matrix) for maximazing flux-carrying high genes
MILPproblem = addGeneMaximization(MILPproblem, HGeneList,RHNames,RLNames, model);
% then, make the gene-centric objective function
% Creating c (objective function)
c_v = zeros(size(S,2),1);
c_yh1 = zeros(length(RHindex),1);
c_yl = ones(length(RLindex),1);
c_yh2 = zeros(length(RHindex),1);
c_gene = ones(length(HGeneList),1);
c = [c_v;c_yh1;c_yl;c_yh2;c_gene];
MILPproblem.c = c;
%% Step 4: solve the MILP and generate the output
solution = solveCobraMILP_XL(MILPproblem, 'timeLimit', runtime, 'logFile', logfile, 'printLevel', verbose);
if solution.stat ~= 1
    error('MILP solving failed! Please inspect the reason!');
end
fprintf('...total high gene fitted: %d \n',sum(solution.int(end-length(c_gene)+1:end)));
yH = boolean(solution.int(1:length(RHindex)) + solution.int((length(RHindex)+length(RLindex)+1):(2*length(RHindex)+length(RLindex))));
nameH = model.rxns(ismember(model.rxns,RHNames));
OpenedHReaction = nameH(yH);
nameL = model.rxns(ismember(model.rxns,RLNames));
yL = boolean(solution.int((length(RHindex)+1):(length(RHindex)+length(RLindex))));
ClosedLReaction = nameL(yL);
OpenGene = HGeneList(boolean(solution.int((2*length(RHindex)+length(RLindex))+1:end)));
end

