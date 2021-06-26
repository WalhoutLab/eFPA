function [MILPproblem] = addAbsFluxVariables(MILPproblem, model)
% add new variables(colunms in A matrix) at the end of the input matrix to
% proxy the absolute value of the v variables (flux of each reaction)
%
% USAGE:
%   [MILPproblem] = addAbsFluxVariables(MILPproblem, model)
%
% INPUT:
%   MILPproblem:    the MILProblem that the new variables will be added to
%                   (could also be a LP)
%   model:          cobra model structure
%
% OUTPUT:
%   MILPproblem:    the mew MILProblem with absolute flux proxy variables added
%
% ..AUTHOR:   Xuhang Li, March 2020
% create the A matrix
A = sparse(size(MILPproblem.A,1)+2*(length(model.rxns)),size(MILPproblem.A,2)+length(model.rxns));
[m,n,s] = find(MILPproblem.A);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end
% add variables 
for i = 1:length(model.rxns)
    A(size(MILPproblem.A,1)+i,i) = -1;
    A(size(MILPproblem.A,1)+i,size(MILPproblem.A,2)+i) = 1;
end

for i = 1:length(model.rxns)
    A(size(MILPproblem.A,1)+length(model.rxns)+i,i) = 1;
    A(size(MILPproblem.A,1)+length(model.rxns)+i,size(MILPproblem.A,2)+i) = 1;
end

% update other fields
lb = [MILPproblem.lb;zeros(length(model.rxns),1)];
ub = [MILPproblem.ub;1000*ones(length(model.rxns),1)];
b = [MILPproblem.b;zeros(2*length(model.rxns),1)];
if size(MILPproblem.csense,1) == 1 %for some model, the csense is a string instead of vector
    csense1(1:(2*length(model.rxns))) = 'G';
    csense = [MILPproblem.csense,csense1];
else %is a vector
    csense1(1:(2*length(model.rxns)),1) = 'G';
    csense = [MILPproblem.csense; csense1];
end
% modify vartype if it is an MILP
if isfield(MILPproblem,'vartype')
    vartype1(1:length(model.rxns),1) = 'C';
    vartype = [MILPproblem.vartype;vartype1];
    MILPproblem.vartype = vartype;
end
% generate function output
MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.csense = csense;

end