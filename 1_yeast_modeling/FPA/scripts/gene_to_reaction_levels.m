function [reaction_levels, reaction_status] = gene_to_reaction_levels( model, genes, levels, f_and, f_or )
% the GPR parser is built by developing the original parser in Machado,
% Daniel, and Markus Herrgård "Systematic evaluation of methods for integration of transcriptomic data into constraint-based models of metabolism." PLoS computational biology 10, no. 4 (2014).
% -------------------------------------------------
% Convert gene expression levels to reaction levels using GPR associations.
% Level is 0 if there is no GPR for the reaction or no measured genes.
%
% INPUTS
%       model - cobra model
%       genes - gene names
%       levels - gene expression levels
%       f_and - function to replace AND
%       f_or - function to replace OR
%
% OUTPUTS
%       reaction_levels - reaction expression levels
%
% Original Author: Daniel Machado, 2013
% Modified by: Xuhang Li, Mar 2020 
% the output now is a string of expression levels for "AND" connected blocks

    reaction_levels = cell(length(model.rxns), 1);
    reaction_status = 1000*ones(length(model.rxns), 1);
    preParsedGrRules = preparseGPR_xl(model.grRules);
    for i = 1:length(model.rxns)
        [level,status] = eval_gpr(preParsedGrRules{i}, genes, levels, f_and, f_or);
        reaction_levels{i} = level;
        reaction_status(i) = status;
    end
end

