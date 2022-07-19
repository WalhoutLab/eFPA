%% About
% in the general predictor of metabolite abundance, we took the max delta
% rFP of all associated internal reactions for a metabolite of interest.
% This scripts prepares the maximum delta rFP of all associated internal
% reactions for a given metabolite

setEnvForAnalysis
%%
paraStr = 'protein_TS_common_originalFPA_originalDist_order100_naiveNetwork';
% protein_TS_common_newFPA_weightedDist_order6_tissueNetwork
% RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork
% RNA_TS_all_newFPA_weightedDist_order6_tissueNetwork
% protein_TS_common_originalFPA_originalDist_order100_naiveNetwork

% load the FPA result of all internal rxns
load(['output/FPA_rxn_',paraStr,'.mat']);
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan and center by median
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_ctd = normalize(relFP,2,'center','median');

% take the corresponding maximum for each metabolite in the model
allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
relFP_nearest = [];
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
rowlabels_ID = regexprep(rowlabels,'_.$','');
S_logical = full(logical(model.S~=0));
for i = 1:length(allCmp_iHumanName)
    myMet = allCmp_iHumanName{i};
    metInd = strcmp(myMet,metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myRxns = intersect(myRxns,targetRxns); % intersect with all internal rxns
    myInd = ismember(rowlabels_ID,myRxns);
    if any(myInd)
         relFP_nearest(i,:) = max(relFP_ctd(myInd,:),[],1);
    else
        relFP_nearest(i,:)= zeros(1,size(relFP_nearest,2));
        fprintf('%s is not predicted (no internal rxn found)\n',allCmp_iHumanName{i});
    end
end
cmpName_nearest = allCmp_iHumanName;

save(['output/pseudoDM_',paraStr,'.mat'],'relFP_nearest','cmpName_nearest')
    