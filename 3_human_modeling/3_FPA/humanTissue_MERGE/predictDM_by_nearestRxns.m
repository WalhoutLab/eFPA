setEnvForAnalysis
%%
load output_tuning/FPA_rxn_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
          
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;
relFP_wtd_ctd = normalize(relFP_wtd,2,'center','median');

%%
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
    myRxns = intersect(myRxns,targetRxns);
    myInd = ismember(rowlabels_ID,myRxns);
    if any(myInd)
         relFP_nearest(i,:) = max(relFP_wtd_ctd(myInd,:),[],1);
    else
        relFP_nearest(i,:)= zeros(1,size(relFP_nearest,2));
        fprintf('%s is not predicted (no internal rxn found)\n',allCmp_iHumanName{i});
    end
end
cmpName_nearest = allCmp_iHumanName;

save('output_tuning/pseudoDM_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat','relFP_nearest','cmpName_nearest')

   
    
    
    
    
    