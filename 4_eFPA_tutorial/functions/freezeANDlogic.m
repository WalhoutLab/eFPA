function expression_logic = freezeANDlogic(expression)
% "Freeze" the levels for "AND" connected blocks. The principle of the GPR
% parser is to recursively merge all non-functional "AND" connected levels
% by the GPR rule. Once a resulted level becomes functionally "AND"
% connected, we will freeze it (to help it be free from getting further
% merged) by adding a text to the begining and end of the number.
%
% INPUTS
%       expression - the GPR rule string that is converted to expression
%       levels
%
% OUTPUTS
%       expression_logic - the functional "AND" freezed GPR string
%
% Author: Xuhang Li, Mar 2020 

NONETYPE = 'NaN';% by default, the undetected genes are NaN and will be ignored
NUMBER = '[0-9\.\-e]+';
MAYBE_NUMBER = [NUMBER '|' NONETYPE];
% freeze the AND logic genes
expression = regexprep(expression,'and','&');
expression = regexprep(expression,'or','|');
numStarts = regexp(expression,MAYBE_NUMBER);
numEnds = regexp(expression,MAYBE_NUMBER,'end');
judgeExpression = expression; % make a string for logic judgement purpose later
for i = 1:length(numStarts)
    judgeExpression(numStarts(i):numEnds(i)) = '1';
end
if numStarts(1) > 1
    expression_logic = expression(1:numStarts(1)-1);
else
    expression_logic = '';
end
for i = 1:length(numStarts)-1
    % replace the gene to be evaluated with "0" and then judge if it will influence the outcome of the logical expression. (if true, it is  a functional "AND" connected gene/block) 
    tmp = [judgeExpression(1:numStarts(i)-1),'0',judgeExpression(numEnds(i)+1:end)];
    if ~eval(tmp) %if this is "AND" logic gene, freeze
        expression_logic = [expression_logic, 'fz_',expression(numStarts(i):numEnds(i)),'_fz',expression(numEnds(i)+1:numStarts(i+1)-1)];
    else
        expression_logic = [expression_logic, expression(numStarts(i):numStarts(i+1)-1)];
    end
end
% specially treat the last number
if numEnds(end) < length(expression)
    tmp = [judgeExpression(1:numStarts(end)-1),'0',judgeExpression(numEnds(end)+1:end)];
    if ~eval(tmp) % if this is "AND" logic gene
        expression_logic = [expression_logic, 'fz_', expression(numStarts(end):numEnds(end)),'_fz',expression(numEnds(end)+1:end)];
    else
        expression_logic = [expression_logic, expression(numStarts(end):end)];
    end
else
    tmp = [judgeExpression(1:numStarts(end)-1),'0'];
    if ~eval(tmp) % if this is "AND" logic gene
        expression_logic = [expression_logic, 'fz_', expression(numStarts(end):numEnds(end)),'_fz'];
    else
        expression_logic = [expression_logic, expression(numStarts(end):end)];
    end
end
expression_logic = regexprep(expression_logic,'&','and'); %change the symbol for evaluation back to text
expression_logic = regexprep(expression_logic,'\|','or'); %change the symbol for evaluation back to text