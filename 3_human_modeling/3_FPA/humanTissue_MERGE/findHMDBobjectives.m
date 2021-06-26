%% the p-value of boxplots should be changed to rank-sum test
%%
setEnvForAnalysis
%% EX transporter and DM for these metabolites and nearest network
% target DM
load('allCmp_iHumanName.mat');
allCmp_iHumanName = unique(allCmp_iHumanName);
targetRxns_DM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);

% target transporters
targetExRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
targetExRxns = intersect(targetExRxns,EXinvolvedRxns);
targetExRxns = [targetExRxns;{'MI1Pt'}];
targetRxns = targetExRxns;
targetRxns_all_transporter = targetRxns;

metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
metInd = ismember(metNames, allCmp_iHumanName);
myRxns = model.rxns(any(model.S(metInd,:),1));
targetRxns_TSP = intersect(targetRxns_all_transporter,myRxns);

% nearest network reactions for these metabolites 
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
S_logical = full(logical(model.S~=0));
% find the target internal reactions 
exludeSys = {'Transport reactions','Exchange/demand reactions'};
upkrxns = model.rxns(any(model.S(ismember(model.mets,{'majorNutr','sideNutr'}),:),1));
intRxns = model.rxns(~ismember(model.subSystems,exludeSys) & ~ismember(model.rxns,upkrxns));
targetRxns_INT = {};
for i = 1:length(allCmp_iHumanName)
    myMet = allCmp_iHumanName{i};
    metInd = strcmp(myMet,metNames);
    myRxns = model.rxns(any(S_logical(metInd,:),1));
    myRxns = intersect(myRxns,intRxns);
    targetRxns_INT = [targetRxns_INT; myRxns];
end
targetRxns_INT = unique(targetRxns_INT);

save('HMDB_objectives.mat','targetRxns_DM','targetRxns_TSP','targetRxns_INT');
   