function [level_f, level_r] = sumExp(targetRxn,n, labels, distMat, valid_rxns, normalizedLevel)

level_f = NaN;
level_r = NaN;

% forward direction
if any(strcmp(labels, [targetRxn,'_f']))
    myDist = distMat(strcmp(labels, [targetRxn,'_f']),:);
    % find the measured reactions within the boundary 
    rxnSet = labels(myDist <= n);
    rxnSet = regexprep(rxnSet,'_.$','');
    % take the average expression 
    if any(ismember(valid_rxns, rxnSet))
        level_f = mean(normalizedLevel(ismember(valid_rxns, rxnSet)));
    else
        level_f = 1;
    end
end
if any(strcmp(labels, [targetRxn,'_r']))
    myDist = distMat(strcmp(labels, [targetRxn,'_r']),:);
    % find the measured reactions within the boundary 
    rxnSet = labels(myDist <= n);
    rxnSet = regexprep(rxnSet,'_.$','');
    % take the average expression 
    if any(ismember(valid_rxns, rxnSet))
        level_r = mean(normalizedLevel(ismember(valid_rxns, rxnSet)));
    else
        level_r = 1;
    end
end
if isnan(level_f)
    level_f = 0;
end
if isnan(level_r)
    level_r = 0;
end

end