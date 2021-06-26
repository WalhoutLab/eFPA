function [result, status] = eval_gpr(rule, genes, levels, f_and, f_or)
% the GPR parser is built by developing the original parser in Machado,
% Daniel, and Markus Herrgård "Systematic evaluation of methods for 
% integration of transcriptomic data into constraint-based models of 
% metabolism." PLoS computational biology 10, no. 4 (2014).
% -------------------------------------------------
% Convert gene expression levels to reaction levels using GPR associations.
% Level is 0 if there is no GPR for the reaction or no measured genes.
%
% INPUTS
%       rule - the GPR rule string
%       genes - gene names
%       levels - gene expression levels
%       f_and - function to replace AND
%       f_or - function to replace OR
%
% OUTPUTS
%       result - a string of expression levels for "AND" connected blocks
%
% Original Author: Daniel Machado, 2013
% Modified by: Xuhang Li, Mar 2020 
EVAL_OK = 1;
PARTIAL_MEASUREMENTS = 0;
NO_GPR_ERROR = -1;
NO_MEASUREMENTS = -2;
MAX_EVALS_EXCEEDED = -3;

MAX_EVALS = 1000;
NONETYPE = 'NaN';% by default, the undetected genes are NaN and will be ignored;

NUMBER = '[0-9\.\-e]+';
MAYBE_NUMBER = [NUMBER '|' NONETYPE];

expression = rule;
result = 'NaN';
status = EVAL_OK;

if isempty(expression)
    status = NO_GPR_ERROR;
else
    % we allow letters, numbers, dot, dash and colon in gene names. Any
    % other special symbol needs to be added below.
    rule_genes = setdiff(regexp(expression,'\<(\w|\-|\.|:)+\>','match'), {'and', 'or'});

    total_measured = 0;
    
    rule_genes = sort(rule_genes); % order the genes from longest to shortest if there are overlapped names. Fix some bugs in rare cases
    rule_genes = rule_genes(length(rule_genes):-1:1);
    
    for i = 1:length(rule_genes)
        j = find(strcmp(rule_genes{i}, genes));
        if isempty(j)
            level = NONETYPE;
        else
            level = num2str(sum(levels(j)));% if a gene appears multiple time, we use the sum of levels!
            total_measured = total_measured + 1;
            last_valid_level = level;
        end
        expression = regexprep(expression, ['\<', rule_genes{i}, '\>'], level );
    end
    if total_measured < length(rule_genes)
        status = PARTIAL_MEASUREMENTS;
    end
    if total_measured > 1 % processing multiple-gene GPR string
        expression_logic = freezeANDlogic(expression);
        maybe_and = @(a,b)maybe_functor(f_and, a, b);
        maybe_or = @(a,b)maybe_functor(f_or, a, b); 
        str_wrapper = @(f, a, b)num2str(f(str2double(a), str2double(b)));

        counter = 0;

        % fold all the "OR"-connected genes
        while contains(expression_logic,'or') 
            counter = counter + 1;
            if counter > MAX_EVALS
                error('parsing failed for grRule "%s".\nPlease check the reason!',rule)
                break
            end
            paren_expr = ['\(\s*(', MAYBE_NUMBER,')\s*\)'];
            and_expr = ['(',MAYBE_NUMBER,')\s+and\s+(',MAYBE_NUMBER,')'];
            or_expr = ['(',MAYBE_NUMBER,')\s+or\s+(',MAYBE_NUMBER,')'];
            len_pre = length(expression_logic);
            expression_logic = regexprep(expression_logic, and_expr, '${str_wrapper(maybe_and, $1, $2)}');
            expression_logic = regexprep(expression_logic, or_expr, '${str_wrapper(maybe_or, $1, $2)}');
            len_post = length(expression_logic);
            if len_pre == len_post %wrap needed
                expression_logic = regexprep(expression_logic, paren_expr, '$1');
            end
            expression_logic = freezeANDlogic(regexprep(expression_logic,'[fz_|_fz]','')); % freeze the unfolded "AND"
        end
        result = regexprep(expression_logic,'[fz_|_fz|(|)]','');
        result = regexprep(result,' +',' ');

    elseif total_measured == 0
        status = NO_MEASUREMENTS;
        result = 'NaN';
    else % only one valid measurement; just put the measurement there
        % all GPR is ignored!
        result = last_valid_level;
    end
end

% post processing ==> convert to the matrix in cell
result = str2double(strsplit(result,' and '));
end
function c = maybe_functor(f, a, b)
    
    if isnan(a) && isnan(b)
        c = nan;
    elseif ~isnan(a) && isnan(b)
        c = a;
    elseif isnan(a) && ~isnan(b)
        c = b;
    else 
        c = f(a,b);
    end
end
