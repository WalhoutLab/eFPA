function [MILProblem] = solution2constraint(MILProblem,solution)
% convert the objective value in the provided solution to a constraint in
% the input MILProblem. Work with both integer objective value or continous
% objective value. Note that the osense value and c vector in the input MILP is assumed
% to be the same as those for the provided solution. 
%
% USAGE:
%   [MILProblem] = solution2constraint(MILProblem,solution)
%
% INPUTS:
%   MILProblem:     COBRA MILP struct (to be added constraint on)
%   solution:       the MILP solution struct for providing the objective
%                   value. Only .obj field is required.
%
% OUTPUTS:
%   MILProblem:     the constrianed MILP. A new row is added to the A
%                   matrix that constrians the MILP to meet the objective value in the solution
%
% AUTHOR: Xuhang Li, Mar 2020

obj = solution.obj;
A = MILProblem.A;
% add constraint
A(end+1,:) = MILProblem.c;
MILProblem.A = A;
% update others fields
MILProblem.b = [MILProblem.b;obj];
if MILProblem.osense == -1
    mycsense = 'G';
else
    mycsense = 'L';
end
if size(MILProblem.csense,1) == 1 % for some model, the csense is a string instead of vector
    csense = [MILProblem.csense,mycsense];
else % when it is a vector
    csense = [MILProblem.csense; mycsense];
end
MILProblem.csense = csense;
end