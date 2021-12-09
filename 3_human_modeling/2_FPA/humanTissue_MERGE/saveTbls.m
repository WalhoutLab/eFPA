%% About 
% produce the output tables for human tissue modeling
% we consider the predictions made by protein data as the final prediction.

setEnvForAnalysis
%% load data and write the flux potential predictions of tissue flux activity
load output/FPA_rxn_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
          
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;

load output/FPA_rxn_RNA_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
relFP(rmInd,:) = [];
relFP_RNA = relFP;

% FC cutoff
relFP_filtered = relFP_wtd - median(relFP_wtd,2);
rowlabels_filtered = rowlabels;
relFP_filtered_RNA = relFP_RNA - median(relFP_RNA,2);

% filter by OFD consistency
lib = 'protein';
for i = 1:length(conditions)
    load(['./../../1_iMAT++/output/humanModel/',lib,'/OFD/',conditions{i},'.mat']);
    eval([conditions{i},'=myCSM;']);
end
IsInOFD = false(length(rowlabels_filtered),1);
for j = 1:length(rowlabels_filtered)
    highestSites = conditions((relFP_filtered(j,:) > 0.2));
    for i = 1:length(highestSites)
        eval(['myCSM=',highestSites{i},';']);
        myrxn = rowlabels_filtered{j};
        myrxnID = myrxn(1:end-2);
        if any(strcmp(myrxnID,targetRxns)) % reaction centric
            if strcmp(myrxn(end),'f')
                if myCSM.OFD(strcmp(model.rxns,myrxnID)) > 0
                    IsInOFD(j) = true;
                    break
                end
            else
                if myCSM.OFD(strcmp(model.rxns,myrxnID)) < 0
                    IsInOFD(j) = true;
                    break
                end
            end
        else % metabolite centric
            if any(strcmp(model.rxns,myrxnID)) % existed rxns
                if strcmp(myrxn(end),'f')
                    if levels_f(strcmp(model.rxns,myrxnID)) >= 0
                        IsInOFD(j) = true;
                        break
                    end
                else
                    if levels_r(strcmp(model.rxns,myrxnID)) >= 0
                        IsInOFD(j) = true;
                        break
                    end
                end
            else % newly added metabolite demand; not comparable with OFD
                IsInOFD(j) = true;
                break
            end
        end
    end
end
IsInOFD_pro = IsInOFD;

lib = 'RNA';
for i = 1:length(conditions)
    load(['./../../1_iMAT++/output/humanModel/',lib,'/OFD/',conditions{i},'.mat']);
    eval([conditions{i},'=myCSM;']);
end
% select by high sites in OFD
IsInOFD = false(length(rowlabels_filtered),1);
for j = 1:length(rowlabels_filtered)
    highestSites = conditions((relFP_filtered_RNA(j,:) > 0.2));
    for i = 1:length(highestSites)
        eval(['myCSM=',highestSites{i},';']);
        myrxn = rowlabels_filtered{j};
        myrxnID = myrxn(1:end-2);
        if any(strcmp(myrxnID,targetRxns)) % reaction centric
            if strcmp(myrxn(end),'f')
                if myCSM.OFD(strcmp(model.rxns,myrxnID)) > 0
                    IsInOFD(j) = true;
                    break
                end
            else
                if myCSM.OFD(strcmp(model.rxns,myrxnID)) < 0
                    IsInOFD(j) = true;
                    break
                end
            end
        else % metabolite centric
            if any(strcmp(model.rxns,myrxnID)) % existed rxns
                if strcmp(myrxn(end),'f')
                    if levels_f(strcmp(model.rxns,myrxnID)) >= 0
                        IsInOFD(j) = true;
                        break
                    end
                else
                    if levels_r(strcmp(model.rxns,myrxnID)) >= 0
                        IsInOFD(j) = true;
                        break
                    end
                end
            else % newly added metabolite demand; not comparable with OFD
                IsInOFD(j) = true;
                break
            end
        end
    end
end
IsInOFD_RNA = IsInOFD;

% label out the high-confident sites in OFD (based on protein data only)
for i = 1:length(rowlabels_filtered)
    if IsInOFD_pro(i)
        rowlabels_filtered{i} = [rowlabels_filtered{i},'*'];
    end
end

% save the table
IDs = regexprep(conditions,'_','-');
relFluxPotential = array2table(relFP_filtered);
relFluxPotential.Properties.VariableNames = IDs;
relFluxPotential.Properties.RowNames = rowlabels_filtered;
writetable(relFluxPotential,'output/supp_delta_rFP_rxn_protein.csv','WriteRowNames',true);

%% laod data and save the delta rFP predictions (general predictor) for tissue-enriched metabolites
clear
%%
setEnvForAnalysis

% define reaction sets
% the exchange reactions with environment 
excRxns = model.rxns(findExcRxns_XL(model));
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
excRxns = intersect(excRxns,EXinvolvedRxns);

% the transporter (with env) reactions
allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
TSP = [];
for i = 1:length(allCmp_iHumanName)
    myMet_e = {[allCmp_iHumanName{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,allCmp_iHumanName(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    % we define the transporter as the reactions that contain the same
    % metabolite in [e] and another compartment (cellular) in the same
    % reaction
    candidate = intersect(myRxns_non_e,myRxns_e);
    % check if is on diff side of the reaction
    if ~isempty(candidate)
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                TSP = union(TSP,candidate(j));
            end
        end
    end
end
% some special transporter will be missed, we add back 
envTspRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
envTspRxns = intersect(envTspRxns,EXinvolvedRxns);
tspRxns = union(TSP, envTspRxns);

% the internal transporters 
% internal transporters are hard to define, since a reaction can span two
% compartments internally but it is not a real transporter. So, we use the
% subsys annotation as a compromise
intTspRxns = setdiff(model.rxns(ismember(model.subSystems,{'Transport reactions'})),tspRxns);

% internal regular reactions
intRxns = setdiff(model.rxns, [excRxns; tspRxns; intTspRxns]);

% regular met-analysis target (transporters)
targetRxns = tspRxns;

load('allCmp_iHumanName.mat');
allCmp_iHumanName = unique(allCmp_iHumanName);
targetRxns_DM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);

allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
targetRxns_allMetDM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);

% load datasets 
% protein prediction
load output/FPA_transporter_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];   
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;
load output/FPA_demand_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_DM = [cellfun(@(x) [x,'_f'],targetRxns_DM,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns_DM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_DM(rmInd) = [];
relFP_wtd_DM = relFP;
load output/pseudoDM_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
cmpName_nearest_network = cmpName_nearest;
relFP_nearest_network = relFP_nearest;
load output/FPA_allMetDemand_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM(rmInd) = [];
relFP_wtd_allMetDM = relFP;

% delta rFP matrix - general predictor (transporter or DM or network-FPA-nearest) - protein 
relFP_sel_allMetDM = zeros(size(relFP_wtd_allMetDM,1),length(measuredTissue));
relFP_wtd_ctd_allMetDM = normalize(relFP_wtd_allMetDM,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_allMetDM(:,i) = max(relFP_wtd_ctd_allMetDM(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
predMat_all = [];
qryMets = unique(metNames);
rowlabels_ID = regexprep(rowlabels,'_.$','');
S_logical = full(logical(model.S~=0));
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);
    
    DMrxn = {['NewMet_',qryMets{j},'_f'],['NewMet_',qryMets{j},'_r']};
    DMInd = ismember(rowlabels_allMetDM,DMrxn);
    
    NearestInd = ismember(cmpName_nearest_network,qryMets{j});
    
    if any(myInd) || any(DMInd) || any(NearestInd)
        predMat_all = [predMat_all;max([relFP_sel(myInd,:);relFP_sel_allMetDM(DMInd,:);relFP_sel_nearest(NearestInd,:)],[],1)];
    else
        predMat_all = [predMat_all;nan(1,size(predMat_all,2))];
    end
end

predMat_all(abs(predMat_all)>1.01) = NaN;% remove numeric errors
% save the table
relFluxPotential = array2table(predMat_all);
relFluxPotential.Properties.VariableNames = measuredTissue;
relFluxPotential.Properties.RowNames = qryMets;
writetable(relFluxPotential,'output/supp_delta_rFP_met_protein.csv','WriteRowNames',true);



