%% About
% convert the rFP predictions to the primary and secondary
% tissue-enrichment sites 
%% set up the env variables
addpath ~/cobratoolbox/
addpath ./input/
addpath ./../../scripts/
addpath ./../../scripts/oriMERGE/
addpath ./../../../input/GEMs/
addpath ./../../../bins/
load('Tissue.mat');
expressionTbl = readtable('expressionTable.tsv','FileType','text','ReadRowNames',true);
tissueLabel = expressionTbl.Properties.VariableNames;
% In this script, we reproduce the reactions in table S7 and S8
FPAtbl = readtable('relativeFluxPotentials.tsv','FileType','text','ReadRowNames',true);% same as table S7 and S8
alltarget = unique(cellfun(@(x) x(1:end-1),FPAtbl.Properties.RowNames,'UniformOutput',false));
targetExRxns = alltarget(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE|SNK|DMN)','once')),alltarget));
targetMets = alltarget(cellfun(@(x) ~isempty(regexp(x,'\[.\]','once')),alltarget)); % targets for metabolite-centric calculation
targetMets2 = cellfun(@(x) ['NewMet_',x],targetMets,'UniformOutput',false);
% merge Ex rxns with Met (give a marker)
targetExRxns = [targetExRxns;targetMets2];
targetRxns = setdiff(alltarget,union(targetMets,targetExRxns));
targetRxns = [targetRxns; alltarget(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE|SNK|DMN)','once')),alltarget));targetMets];

%% write out result
algorithms = {'oriMERGE','wtdDist','wtdDist_exp_decay_base2'};
for i = 1:length(algorithms)
    load(['output/FPA_rxns_',algorithms{i},'.mat']);
    rFP = [relFP_f;relFP_r];
    labels = [cellfun(@(x) [x,'_f'], targetRxns,'UniformOutput',0);cellfun(@(x) [x,'_r'], targetRxns,'UniformOutput',0)];
    % remove not applicable
    ind1 = all(isnan(rFP),2);
    rFP(ind1,:) = [];
    labels(ind1) = [];
    % put asside problematic calculations if any
    ind1 = any(isnan(rFP),2);
    rxns_error = labels(ind1);
    rFP(ind1,:) = [];
    labels(ind1) = [];
    % minmax filter
    ind1 = (max(rFP,[],2) - min(rFP,[],2)) >= 0.35;
    rFP = rFP(ind1,:);
    labels = labels(ind1);
    % sorting algorithm
    primarySites = repmat({'NA'} , length(labels),1);
    SecondSites = repmat({'NA'} , length(labels),1);
    for j = 1:size(rFP,1)
        FPvect = rFP(j,:);
        tissues = tissueLabel;
        [FPvect order] = sort(FPvect,'descend');
        tissues = tissues(order);
        if (FPvect(1) >= 1.4*median(FPvect) || FPvect(1) >= median(FPvect) + 0.2) && ...
                FPvect(1) >= FPvect(3) + 0.1
            primarySites(j) = tissues(1);
        end
        if (FPvect(2) >= 1.2*median(FPvect) || FPvect(2) >= median(FPvect) + 0.1) && ...
                FPvect(2) >= FPvect(3) + 0.07
            SecondSites(j) = tissues(2);
        end
    end
    save(['output/FPA_interpretation_',algorithms{i},'.mat'],'labels','primarySites','SecondSites','rxns_error');
end
