function [Flux_min]= minimizeModelFlux_XL(model)
% THIS FUNCTION IS NOT USED IN STANDARD MERGE PIPELINE!
% FOR DEBUGING USE ONLY
% This is a flux minimization tool that will minimize the total flux of a
% input model, regardless of the model format (so the csense could be
% flexible, i.e, contains unequal ones)
MILPproblem.A = model.S;
MILPproblem.lb = model.lb;
MILPproblem.ub = model.ub;
MILPproblem.b = model.b;
MILPproblem.csense = model.csense;
if isfield(model,'vartype')
    MILPproblem.vartype = model.vartype;
end
MILPproblem.c = model.c;

A = sparse(size(MILPproblem.A,1)+2*(length(model.rxns)),size(MILPproblem.A,2)+length(model.rxns));
[m,n,s] = find(MILPproblem.A);
for i = 1:length(m)
    A(m(i),n(i)) = s(i);
end
% add flux measuremnet 
for i = 1:length(model.rxns)
    A(size(MILPproblem.A,1)+i,i) = -1;
    A(size(MILPproblem.A,1)+i,size(MILPproblem.A,2)+i) = 1;
end
for i = 1:length(model.rxns)
    A(size(MILPproblem.A,1)+length(model.rxns)+i,i) = 1;
    A(size(MILPproblem.A,1)+length(model.rxns)+i,size(MILPproblem.A,2)+i) = 1;
end
MILPproblem.A = A;
%create other inputs
MILPproblem.lb = [MILPproblem.lb;zeros(length(model.rxns),1)];
MILPproblem.ub = [MILPproblem.ub;2000*ones(length(model.rxns),1)];
MILPproblem.b = [MILPproblem.b;zeros(2*length(model.rxns),1)];
if size(MILPproblem.csense,1) == 1 %for some model, the csense is a string instead of vector
    csense1(1:2*length(model.rxns)) = 'G';
    csense = [MILPproblem.csense,csense1];
else %is a vector
    csense1(1:(2*length(model.rxns)),1) = 'G';
    csense = [MILPproblem.csense; csense1];
end
MILPproblem.csense = csense;
% Creating c (objective function)
c_old = zeros(length(MILPproblem.c),1);
c_minFlux = ones(length(model.rxns),1);
c = [c_old;c_minFlux];
MILPproblem.c = c;
% set the objective sense as minimize 
MILPproblem.osense = 1;
% sovle it
solution = solveCobraLP(MILPproblem);
if solution.stat == 1
    Flux_min = solution.full(1:length(model.rxns));
else
    Flux_min = solution;
end
end
