function blockList = getBlockList(model,levelTbl_f,levelTbl_r)
% get the rxn IDs to block in FPA from the FVA result in IMAT++ analysis
%
% USAGE:
%   blockList = getBlockList(model,levelTbl_f,levelTbl_r)
%
% INPUT:
%   model:       cobra model structure
%   levelTbl_f:  the parsed FVA levels of all reactions (forward direction) 
%   levelTbl_r:  the parsed FVA levels of all reactions (reverse direction) 
%
% OUTPUT:
%   blockList:   reaction names to block in FPA
%
% AUTHORS:  Xuhang Li, March 2020
model_irrev = convertToIrreversible(model);
% unify nomenclature
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); % has no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
for i = 1:size(levelTbl_f,2)
    blockList{i} = {};
    for j = 1: length(model_irrev.rxns)
        rxnID = model_irrev.rxns{j};
        rxnID_ori = rxnID(1:end-2);
        direction = rxnID(end);
        if strcmp(direction,'f')
            if levelTbl_f(strcmp(model.rxns,rxnID_ori),i) == -1
                blockList{i} = [blockList{i},{rxnID}];
            end
        else
            if levelTbl_r(strcmp(model.rxns,rxnID_ori),i) == -1
                blockList{i} = [blockList{i},{rxnID}];
            end
        end
    end
end
blockList{end+1} = {};% for supercond (super condition is not blocked)
end
