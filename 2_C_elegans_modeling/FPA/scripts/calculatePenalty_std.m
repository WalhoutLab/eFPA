function [normalizedLevel, normalizedLevel_std] = calculatePenalty_std(model,master_expression,master_expression_std)
% calculate the penalty matrix given expression data. It translates the
% expression level of genes to level of reactions by a set of GPR parsing
% rules. In brief, we taking the "sum" of all "OR" gated genes in a single
% GPR and "min" of "AND" ones. Complex GPRs are handled specially by
% converting the GPR annotation to functionally "AND" gated blocks, where
% all "OR" gated genes are added up. In the convertion to blocks, if a
% "AND" gate is nested in an "OR", the minimal value will be taken. For
% detailed explanation, see supplemental information. 
% 
% USAGE:
%
%    penalty = calculatePenalty(model,master_expression,manualPenalty)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    master_expression: the expression profiles of queried conditions. The
%                       master expression variable is reqired to be a cell
%                       array of structure variables. Each structure
%                       variable corresponds to the expression profile of a
%                       condition in comparison. The structure variable
%                       should have two fields, "genes" and "value"
%                       respectively. For instructions on forming the
%                       master expression variable from a expression table
%                       (i.e, TPM table in text format), please see
%                       "FPA_walkthrough_generic.m"
% OPTIONAL INPUTS:
%    manualPenalty:     the user defined penalty for specific reactions.
%                       This input will overide all penalty calculation based on expression
%                       data.
%
%
% OUTPUT:
%   penalty:            the full penalty matrix to be directly used in FPA
%                       calculation
%
% `Yilmaz et al. (2020). Final Tittle and journal.
%
% .. Author: - (COBRA implementation) Xuhang Li, Mar 2020
%% step1: mapping the gene-centric expression levels to reactions
levels = cell(length(model.rxns),length(master_expression));
status = -1*ones(length(model.rxns),length(master_expression));
fprintf('mapping expression levels to penalties...\n');
fprintf(['\n' repmat('.',1,length(master_expression)) '\n\n']);
parfor i = 1:length(master_expression)
    expression = master_expression{i};
    % expression.value = expression.value + 1; %add 1 pseudocount to avoid numerical error
    [levels(:,i), status(:,i)] = gene_to_reaction_levels(model, expression.genes, expression.value, @min, @(x,y)(x+y));%GPR parser to convert the expression of genes to levels of "AND" gated blocks
    fprintf('\b|\n');%for simple progress monitor
end

levels_std = cell(length(model.rxns),length(master_expression_std));
status = -1*ones(length(model.rxns),length(master_expression_std));
fprintf('mapping expression levels to penalties...\n');
fprintf(['\n' repmat('.',1,length(master_expression_std)) '\n\n']);
parfor i = 1:length(master_expression_std)
    expression = master_expression_std{i};
    % expression.value = expression.value + 1; %add 1 pseudocount to avoid numerical error
    [levels_std(:,i), status(:,i)] = gene_to_reaction_levels(model, expression.genes, expression.value, @min, @(x,y)(x+y));%GPR parser to convert the expression of genes to levels of "AND" gated blocks
    fprintf('\b|\n');%for simple progress monitor
end
%% step2: calculate penalty
normalizedLevel = nan(size(levels,1),size(levels,2));%the levels normalized to super condition
penalty = ones(size(levels,1),size(levels,2));%set an default penalty as 1
normalizedLevel_std = nan(size(levels,1),size(levels,2));%the levels normalized to super condition
penalty_std = ones(size(levels,1),size(levels,2));%set an default penalty as 1
for i = 1:length(model.rxns)
    % this is for error handling purposes. If some very special GPR format
    % appears and the parser fails to correctly handle it, the number of
    % blocks for the same reaction in different conditions is likely
    % different. So, we check the consistance first to make sure the GPR
    % parsing is as expected.
    lenL = zeros(length(master_expression));
    for j=1:length(master_expression)
        lenL(j)=length(levels{i,j});
    end
    if any(lenL ~= lenL(1)) %some GPR parsing error
        error('GPR length unequal');
    else
        % normalize every block
        stackM = nan(lenL(1),length(master_expression)); %stack the levels of blocks in a matrix (it was in cell variable originally)
        stackM_std = nan(lenL(1),length(master_expression)); %stack the levels of blocks in a matrix (it was in cell variable originally)
        for z = 1:lenL(1)
            for s = 1:length(master_expression)
                stackM(z,s) = levels{i,s}(z);
                stackM_std(z,s) = levels_std{i,s}(z);
            end
        end
        % first normalize the matrix to the super condition (highest value)
        stackM_std = stackM_std ./ max(stackM,[],2);
        stackM = stackM ./ max(stackM,[],2);
        % then, pick the minimal value of all blocks (because they are "AND" connected
        [tmp minInd] = min(stackM,[],1);
        normalizedLevel(i,:) = tmp;
        normalizedLevel_std(i,:) = stackM_std(sub2ind(size(stackM_std),minInd,1:size(stackM_std,2)));
        if ~any(isnan(normalizedLevel(i,:))) %if no NaN occurs
            penalty(i,:) = 1./normalizedLevel(i,:);
        % if all the levels (in different conditions) are NaN, this
        % reaction must have no data. We will keep the default penalty (which is 1) for
        % it 
        elseif any(isnan(normalizedLevel(i,:))) && ~all(isnan(normalizedLevel(i,:))) %if only some condition gives NaN
            error('partial None Penalty for some conditions, check!');
        end
    end
end

end