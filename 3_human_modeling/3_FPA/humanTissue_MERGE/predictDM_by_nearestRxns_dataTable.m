setEnvForAnalysis
%%
load output/FPA_rxn_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
          
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;
relFP_wtd_ctd = normalize(relFP_wtd,2,'center','median');

logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
% align the tables
colLabels = logTPMTbl.Properties.VariableNames(2:end);
%%
load('allCmp_iHumanName.mat');
allCmp_iHumanName = unique(allCmp_iHumanName);
relFP_nearest = [];
relFP_nearest_detail = {};
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
        detail = array2table(relFP_wtd_ctd(myInd,:));
        detail.Properties.RowNames = rowlabels(myInd);
        detail.Properties.VariableNames = colLabels;
        for j = 1:length(colLabels)
            relFP_nearest_detail(i,j) = {detail(:,j)};
        end
    else
        relFP_nearest(i,:)= zeros(1,size(relFP_nearest,2));
        relFP_nearest_detail(i,:) = repmat({''},1,length(colLabels));
        fprintf('%s is not predicted (no internal rxn found)\n',allCmp_iHumanName{i});
    end
end
cmpName_nearest = allCmp_iHumanName;

t_delta_relFP_nearest = array2table(relFP_nearest);
t_delta_relFP_nearest.Properties.RowNames = cmpName_nearest;
t_delta_relFP_nearest_detail = cell2table(relFP_nearest_detail);
t_delta_relFP_nearest_detail.Properties.RowNames = cmpName_nearest;

t_delta_relFP_nearest.Properties.VariableNames = colLabels;
t_delta_relFP_nearest_detail.Properties.VariableNames = colLabels;
save('output/pseudoDM_detailTbl_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat','t_delta_relFP_nearest','t_delta_relFP_nearest_detail')
 
    
    