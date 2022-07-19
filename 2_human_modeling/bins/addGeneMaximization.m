function [MILPproblem] = addGeneMaximization(MILPproblem, HGeneList, HRxnList,LRxnList, model)
% add the new variables (columns in the A matrix) to a MILProblem for doing gene-centric IMAT++ fitting.
% the binary variables G(i) are added such that:
%       G(i) <= R(i,1) + R(i,2) + ... + R(i,n)
%               where G(i) is the binary variable for Gene i, and R(i,1) to
%               R(i,n) are binary variables proxying the flux state for all
%               the high reactions associated with gene i
%
% USAGE:
%   [MILPproblem] = addGeneMaximization(MILPproblem, HGeneList, HRxnList,LRxnList, model)
%
% INPUT:
%   MILPproblem:    the MILProblem that the new variables will be added to
%   HGeneList:      the list of highly expressed genes
%   HRxnList:       the list of high reactions (associated with highly
%                   expressed genes, but defined by GPR mapping)
%   LRxnList:       the list of rarely expressed reactions 
%   model:          cobra model structure
%
% OUTPUT:
%   MILPproblem:    the mew MILProblem with gene-centric optimization
%                   variables added
%
% WARNING:          we only recommand to use this function inside of
%                   IMAT++ optimization package. Since A matrix is assumed to be in the
%                   format of IMAT++ MILP, errors may occur if applied to other MILP
%                   directly.
%
% ..AUTHOR:   Xuhang Li, March 2020

% find the HRxn index in input matrix 
HRindex0 = find(ismember(model.rxns, HRxnList));
% create A matrix
A = sparse(size(MILPproblem.A,1)+length(HGeneList),size(MILPproblem.A,2)+length(HGeneList));
[m,n,s] = find(MILPproblem.A);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end
% add the gene-centric variables
for i = 1:length(HGeneList)
    myRxns = model.rxns(any(model.rxnGeneMat(:,strcmp(HGeneList(i),model.genes)),2),1);
    HRNames = HRxnList(ismember(HRxnList,myRxns));
    HRindex1 = find(ismember(model.rxns, HRNames));
    HRindex = find(ismember(HRindex0,HRindex1));
    A(size(MILPproblem.A,1)+i,length(model.rxns)+HRindex) = 1;
    A(size(MILPproblem.A,1)+i,length(model.rxns)+length(HRxnList)+length(LRxnList)+HRindex) = 1;
    A(size(MILPproblem.A,1)+i,size(MILPproblem.A,2)+i) = -1;
end

% update other MILP fields
lb = [MILPproblem.lb;zeros(length(HGeneList),1)];
ub = [MILPproblem.ub;ones(length(HGeneList),1)];
b = [MILPproblem.b;zeros(length(HGeneList),1)];
csense1(1:length(HGeneList)) = 'G';
csense = [MILPproblem.csense,csense1];
vartype1(1:length(HGeneList),1) = 'B';
vartype = [MILPproblem.vartype;vartype1];

% generate output MILP
MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.csense = csense;
MILPproblem.vartype = vartype;
end