%% About
% compare the transporting flux potential of metabolites with no
% integration (target transporter expression only) in their predictive
% power of tissue-enriched metabolites
%% load environment
addpath ./../scripts/
addpath input/
setEnvForAnalysis
addpath('PlotPub/lib')

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
%% load datasets 
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

% control - no network
load output/FPA_transporter_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_noNetwork = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_noNetwork(rmInd) = [];
relFP_wtd_noNetwork = relFP;
load output/pseudoDM_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
cmpName_nearest_noNetwork = cmpName_nearest;
relFP_nearest_noNetwork = relFP_nearest;

%% select the metabolite associated with a transporter measurements 
nochange = all(abs(relFP_wtd_noNetwork - relFP_wtd_noNetwork(:,1))<1e-6,2);
validTspNames = rowlabels_noNetwork(~nochange);
validTspNames = regexprep(validTspNames,'_(f|r)$','');
% ==> at least partial expression information is available for 2027/3663
% transporters
validTspMets = model.metNames(any(model.S(:,ismember(model.rxns,validTspNames))~=0,2));
validTspMets = unique(regexprep(validTspMets,' \[.+\]$',''));
%% load HMDB reference set 
metTissueTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','met_tissue_set_processed');
TissueAligTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','tissueAlignment');
cmpTbl = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','metName');
allCmp = {};
validTissues = {};
for i = 1:length(TissueAligTbl.HMDBtissues)
    validTissues = [validTissues; strsplit(TissueAligTbl.HMDBtissues{i},'; ')'];
end
validTissues = unique(validTissues);
for i = 1:size(metTissueTbl,1)
    if any(strcmp(metTissueTbl.name{i},validTissues))
        cmplist = metTissueTbl{i,5};
        cmplist = strsplit(cmplist{:},'; ')';
        allCmp = union(allCmp,cmplist);
    end
end
load('metTbl.mat');
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
% match HMDB cmp to iHuman cmp
[A B] = ismember(allCmp,cmpTbl.name); 
allCmp_MSEAname = allCmp(A);
allCmp_HMDB = cmpTbl.hmdb_id(B(A));
[A B] = ismember(allCmp_HMDB,metTbl.HMDB);
allCmp_iHumanName = metTbl.name(B(A));
allCmp_MSEAname = allCmp_MSEAname(A);
specialMatch = readtable('input/MSEA_dataset/met_tissue_set_processed.xlsx','Sheet','directMatch');
allCmp_iHumanName = [allCmp_iHumanName;specialMatch.ihumanName];
allCmp_MSEAname = [allCmp_MSEAname;specialMatch.metName];

% we only look at metabolites that have a transporter (the metabolite is being
% transported) between cellular space and extracellular space
% we dont consider the special transporters for this analysis (when the
% cargo is converted during transportation)
isTSP = [];
for i = 1:length(allCmp_iHumanName)
    myMet_e = {[allCmp_iHumanName{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,allCmp_iHumanName(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    candidate = intersect(myRxns_non_e,myRxns_e);
    if ~isempty(candidate)
        pass = 0;
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                pass = 1;
            end
        end
        isTSP(i) = pass;
    else
        isTSP(i) = 0;
    end

end
fprintf('%f.2 of metabolites is transportable\n',sum(isTSP)/length(isTSP));
% additionally, only the transporters with expression data (protein)
isexp = ismember(allCmp_iHumanName,validTspMets)';

allCmp_iHumanName = allCmp_iHumanName(logical(isTSP)&isexp);
allCmp_MSEAname = allCmp_MSEAname(logical(isTSP)&isexp);

% tissue specificity matrix
measuredTissue = TissueAligTbl.name;
refMat = zeros(length(allCmp_iHumanName),length(measuredTissue));
for i = 1:length(measuredTissue)
    metList = {};
    allMets = metTissueTbl.member(ismember(metTissueTbl.name,strsplit(TissueAligTbl.HMDBtissues{i},'; ')));
    for j = 1:length(allMets)
        metList = union(metList,strsplit(allMets{j},'; '));
    end
    refMat(ismember(allCmp_MSEAname,metList),i) = 1;
end
%% delta rFP matrix - protein common - transporter FPA
% first merge the tissues
relFP_sel = zeros(size(relFP_wtd,1),length(measuredTissue));
relFP_wtd_ctd = normalize(relFP_wtd,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel(:,i) = max(relFP_wtd_ctd(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
% then, take max of all associated transporters
predMat = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels,'_.$',''),myRxns);

    if any(myInd)
        predMat(j,:) = max(relFP_sel(myInd,:),[],1);
    else
        predMat(j,:) = zeros(1,size(predMat,2));
        warning('%s is not predicted\n',allCmp_iHumanName{j});
    end
end

%% delta rFP matrix - no network integration - transporter
% same as above
relFP_sel_noNetwork_tsp = zeros(size(relFP_wtd_noNetwork,1),length(measuredTissue));
relFP_wtd_ctd_noNetwork_tsp = normalize(relFP_wtd_noNetwork,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_noNetwork_tsp(:,i) = max(relFP_wtd_ctd_noNetwork_tsp(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
predMat_noNetwork_tsp = [];
for j = 1:length(allCmp_iHumanName)
    metInd = strcmp(allCmp_iHumanName{j},metNames);
    myRxns = model.rxns(any(model.S(metInd,:),1));
    myInd = ismember(regexprep(rowlabels_noNetwork,'_.$',''),myRxns);

    if any(myInd)
        predMat_noNetwork_tsp(j,:) = max(relFP_sel_noNetwork_tsp(myInd,:),[],1);
    else
        predMat_noNetwork_tsp(j,:) = zeros(1,size(predMat_noNetwork_tsp,2));
        warning('%s is not predicted\n',allCmp_iHumanName{j});
    end
end
%% assessing the result: hypergeometric enrichment 
cutoffs = -1:0.01:1;
p_mat =[];
Overall_p =[];
Ncall = [];

p_mat_noNetwork_tsp =[];
Overall_p_noNetwork_tsp =[];
Ncall_noNetwork_tsp = [];
for i = 1: length(cutoffs)
    predictions = predMat >= cutoffs(i) & predMat < 1.01; % filter out numerically problematic predictions
    % test enrichment
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall(i) = sum(N);
    
    predictions = predMat_noNetwork_tsp >= cutoffs(i) & predMat_noNetwork_tsp < 1.01;
    x = sum(predictions & refMat,1);
    M = size(predictions,1);
    K = sum(refMat,1);
    N = sum(predictions,1);
    p_mat_noNetwork_tsp(i,:) = hygecdf(x-1,M,K,N,'upper');
    Overall_p_noNetwork_tsp(i) = hygecdf(sum(x)-1,M .* size(predictions,2),sum(K),sum(N),'upper');
    Ncall_noNetwork_tsp(i) = sum(N);
    
end

figure(1)
hold on
plot(cutoffs, -log10(Overall_p),'-');
plot(cutoffs, -log10(Overall_p_noNetwork_tsp),'-');
hold off
legend({'FPA - transporter','expression - transporter '});
xlabel('cutoff')
ylabel('-log10(P value)')
xlim([-1 1]);
plt = Plot(); % create a Plot object and grab the current figure
plt.LegendLoc = 'NorthWest';
plt.BoxDim = [5.5, 4.5];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export('figures/NoTrack_benchmarkIntegrationBenefit_transporterObj_MeasuredTspOnly.pdf');

%% directly assess the delta rFP
% enrichment of high rFP in the tissue-specific set 
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
qryMets = unique(metNames);% all mets
% here we only look at the expression measured and transportable
% metabolites
isTSP = [];
for i = 1:length(qryMets)
    myMet_e = {[qryMets{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,qryMets(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    candidate = intersect(myRxns_non_e,myRxns_e);
    if ~isempty(candidate)
        pass = 0;
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                pass = 1;
            end
        end
        isTSP(i) = pass;
    else
        isTSP(i) = 0;
    end

end
% additionally, only the transporters with expression data (protein)
isexp = ismember(qryMets,validTspMets)';

fprintf('%f.2 of metabolites is transportable&expressed\n',sum(logical(isTSP)&isexp)/length(isTSP));

qryMets = qryMets(logical(isTSP)&isexp);

%% load background delta rFP matrix - transporter FPA - protein
relFP_sel_allMetDM = zeros(size(relFP_wtd_allMetDM,1),length(measuredTissue));
relFP_wtd_ctd_allMetDM = normalize(relFP_wtd_allMetDM,2,'center','median');
for i = 1:length(measuredTissue)
    relFP_sel_allMetDM(:,i) = max(relFP_wtd_ctd_allMetDM(:,ismember(conditions,strsplit(TissueAligTbl.FPAtissues{i},'; '))),[],2);
end
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
predMat_all = [];
rowlabels_ID = regexprep(rowlabels,'_.$','');
S_logical = full(logical(model.S~=0));
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);

    if any(myInd)
        predMat_all = [predMat_all;max(relFP_sel(myInd,:),[],1)];
    end
end
%% load delta rFP matrix - transporter - no network
predMat_all_noNetwork_tsp = [];
rowlabels_ID = regexprep(rowlabels_noNetwork,'_.$','');
for j = 1:length(qryMets)
    metInd = strcmp(qryMets{j},metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myInd = ismember(rowlabels_ID,myRxns);
    
    DMInd = ismember(cmpName_nearest,qryMets{j});
    if any(myInd) 
        predMat_all_noNetwork_tsp = [predMat_all_noNetwork_tsp;max([relFP_sel_noNetwork_tsp(myInd,:)],[],1)];
    end
end

%% plot the box plots - benchmark integration benefit
p = [];
boxData = [];
boxLabel = [];
xLabels = {};

enrichedMet_rFP_noNetwork_tsp = predMat_noNetwork_tsp(logical(refMat));
boxData = [boxData;enrichedMet_rFP_noNetwork_tsp;predMat_all_noNetwork_tsp(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enrich_noNetwork_tsp'},length(enrichedMet_rFP_noNetwork_tsp),1);repmat({'all_noNetwork_tsp'},size(predMat_all_noNetwork_tsp,1)*size(predMat_all_noNetwork_tsp,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{tranporter expression}';'all metabolites^{tranporter expression}'}];
p(1) = ranksum(enrichedMet_rFP_noNetwork_tsp,predMat_all_noNetwork_tsp(:),  'Tail','right');

enrichedMet_rFP = predMat(logical(refMat));
boxData = [boxData;enrichedMet_rFP;predMat_all(:)];
boxData(abs(boxData)>1.01) = NaN;% invalid due to numeric error
boxLabel = [boxLabel;[repmat({'enriched'},length(enrichedMet_rFP),1);repmat({'all'},size(predMat_all,1)*size(predMat_all,2),1)]];
xLabels = [xLabels;{'enriched metabolites^{FPA}';'all metabolites^{FPA}'}];
p(2) = ranksum(enrichedMet_rFP,predMat_all(:),  'Tail','right');

figure('units','inch','position',[0,0,9,8])
boxplot(boxData, boxLabel,'Labels',xLabels,'Colors','bkbk','Symbol','ro');
ylim([-1.1 1.1])
ylabel('\DeltarFP');
xlabel('');
set(gca, 'TickLabelInterpreter', 'tex');
text(1.3,0.75,['p < ',num2str(p(1),2)],'fontsize', 12,'BackgroundColor','w')
text(3.3,0.75,['p < ',num2str(p(2),2)],'fontsize', 12,'BackgroundColor','w')
set(gca, 'fontsize', 14);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
xtickangle(30)
S = hgexport('readstyle','default_sci');
style.Format = 'pdf';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
saveas(gca,'figures/benchmarkIntegrationBenefit_boxplot_transporterObj_MeasuredTspOnly.pdf');
