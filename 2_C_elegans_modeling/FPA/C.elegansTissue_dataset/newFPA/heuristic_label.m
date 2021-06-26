%% start with some manual inspection
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
algorithms = {'oriMERGE','wtdDist','oriDist_exp_decay_base2','oriDist_exp_decay_base100','wtdDist_exp_decay_base2','wtdDist_exp_decay_base100'};
for i = 1:length(algorithms)
    load(['FPA_rxns_',algorithms{i},'.mat']);
    rFP = [relFP_f;relFP_r];
    labels = [cellfun(@(x) [x,'_f'], targetRxns,'UniformOutput',0);cellfun(@(x) [x,'_r'], targetRxns,'UniformOutput',0)];
    % remove not applicable
    ind1 = all(isnan(rFP),2);
    rFP(ind1,:) = [];
    labels(ind1) = [];
    % put asside problematic calculations 
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
    save(['FPA_interpretation_',algorithms{i},'.mat'],'labels','primarySites','SecondSites','rxns_error');
end

%% QC of the heuristic algorithm: compare with table EV7 EV8
load FPA_interpretation_oriMERGE.mat
EV7 = readtable('Table_EV7.xlsx');
EV7_rxns = EV7.ID;
EV7_rxns = regexprep(EV7_rxns,'f$','_f');
EV7_rxns = regexprep(EV7_rxns,'r$','_r');
primarySites_EV7 = regexprep(EV7.Primary,'\*$','');
primarySites_EV7 = regexprep(primarySites_EV7,'gli','Glia');
primarySites_EV7 = regexprep(primarySites_EV7,'gon','Gonad');
primarySites_EV7 = regexprep(primarySites_EV7,'hyp','Hypodermis');
primarySites_EV7 = regexprep(primarySites_EV7,'int','Intestine');
primarySites_EV7 = regexprep(primarySites_EV7,'mus','Body_wall_muscle');
primarySites_EV7 = regexprep(primarySites_EV7,'neu','Neurons');
primarySites_EV7 = regexprep(primarySites_EV7,'pha','Pharynx');

SecondSites_EV7 = regexprep(EV7.Secondary,'\*$','');
SecondSites_EV7 = regexprep(SecondSites_EV7,'gli','Glia');
SecondSites_EV7 = regexprep(SecondSites_EV7,'gon','Gonad');
SecondSites_EV7 = regexprep(SecondSites_EV7,'hyp','Hypodermis');
SecondSites_EV7 = regexprep(SecondSites_EV7,'int','Intestine');
SecondSites_EV7 = regexprep(SecondSites_EV7,'mus','Body_wall_muscle');
SecondSites_EV7 = regexprep(SecondSites_EV7,'neu','Neurons');
SecondSites_EV7 = regexprep(SecondSites_EV7,'pha','Pharynx');

primarySites_XL = repmat({'NA'} , length(EV7_rxns),1);
SecondSites_XL = repmat({'NA'} , length(EV7_rxns),1);
[A B] = ismember(EV7_rxns,labels);
primarySites_XL(A) = primarySites(B(A));
SecondSites_XL(A) = SecondSites(B(A));

mismatch_pri = [];
mismatch_sec = [];
for i = 1:length(primarySites_EV7)
    if ~strcmp(primarySites_EV7{i},primarySites_XL{i})
        mismatch_pri = [mismatch_pri,i];
    end
    if ~strcmp(SecondSites_EV7{i},SecondSites_XL{i})
        mismatch_sec = [mismatch_sec,i];
    end
end

primarySites_XL(mismatch_pri)
primarySites_EV7(mismatch_pri)

SecondSites_XL(mismatch_sec)
SecondSites_EV7(mismatch_sec)
%% QC of the heuristic algorithm: validate selected cases
load(['FPA_rxns_oriMERGE.mat']);
rFP = [relFP_f;relFP_r];
labels = [cellfun(@(x) [x,'_f'], targetRxns,'UniformOutput',0);cellfun(@(x) [x,'_r'], targetRxns,'UniformOutput',0)];

i = 9;
myrxn = EV7_rxns(mismatch_pri(i));
myrxn2 = regexprep(myrxn,'_','');
FP_EV7 = cellfun(@str2num,([EV7{strcmp(EV7.ID,myrxn2),[8,7,10,13,9,12,11]}]));
FP_XL = rFP(strcmp(labels,myrxn),:);
plot(FP_XL,FP_EV7,'r.')

%% secondary
i =12;
myrxn = EV7_rxns(mismatch_sec(i));
myrxn2 = regexprep(myrxn,'_','');
FP_EV7 = cellfun(@str2num,([EV7{strcmp(EV7.ID,myrxn2),[8,7,10,13,9,12,11]}]));
FP_XL = rFP(strcmp(labels,myrxn),:);
plot(FP_XL,FP_EV7,'r.')

%% EV8
load FPA_interpretation_oriMERGE.mat
EV8 = readtable('Table_EV8.xlsx');
EV8_rxns = EV8.ID;
EV8_rxns = regexprep(EV8_rxns,'f$','_f');
EV8_rxns = regexprep(EV8_rxns,'r$','_r');
EV8_rxns = regexprep(EV8_rxns,']$',']_f');

primarySites_EV8 = regexprep(EV8.Primary,'\*$','');
primarySites_EV8 = regexprep(primarySites_EV8,'gli','Glia');
primarySites_EV8 = regexprep(primarySites_EV8,'gon','Gonad');
primarySites_EV8 = regexprep(primarySites_EV8,'hyp','Hypodermis');
primarySites_EV8 = regexprep(primarySites_EV8,'int','Intestine');
primarySites_EV8 = regexprep(primarySites_EV8,'mus','Body_wall_muscle');
primarySites_EV8 = regexprep(primarySites_EV8,'neu','Neurons');
primarySites_EV8 = regexprep(primarySites_EV8,'pha','Pharynx');

SecondSites_EV8 = regexprep(EV8.Secondary,'\*$','');
SecondSites_EV8 = regexprep(SecondSites_EV8,'gli','Glia');
SecondSites_EV8 = regexprep(SecondSites_EV8,'gon','Gonad');
SecondSites_EV8 = regexprep(SecondSites_EV8,'hyp','Hypodermis');
SecondSites_EV8 = regexprep(SecondSites_EV8,'int','Intestine');
SecondSites_EV8 = regexprep(SecondSites_EV8,'mus','Body_wall_muscle');
SecondSites_EV8 = regexprep(SecondSites_EV8,'neu','Neurons');
SecondSites_EV8 = regexprep(SecondSites_EV8,'pha','Pharynx');

primarySites_XL = repmat({'NA'} , length(EV8_rxns),1);
SecondSites_XL = repmat({'NA'} , length(EV8_rxns),1);
[A B] = ismember(EV8_rxns,labels);
primarySites_XL(A) = primarySites(B(A));
SecondSites_XL(A) = SecondSites(B(A));

mismatch_pri = [];
mismatch_sec = [];
for i = 1:length(primarySites_EV8)
    if ~strcmp(primarySites_EV8{i},primarySites_XL{i})
        mismatch_pri = [mismatch_pri,i];
    end
    if ~strcmp(SecondSites_EV8{i},SecondSites_XL{i})
        mismatch_sec = [mismatch_sec,i];
    end
end

primarySites_XL(mismatch_pri)
primarySites_EV8(mismatch_pri)

SecondSites_XL(mismatch_sec)
SecondSites_EV8(mismatch_sec)

load(['FPA_rxns_oriMERGE.mat']);
rFP = [relFP_f;relFP_r];
labels = [cellfun(@(x) [x,'_f'], targetRxns,'UniformOutput',0);cellfun(@(x) [x,'_r'], targetRxns,'UniformOutput',0)];
%% check mismatches
load(['FPA_rxns_oriMERGE.mat']);
rFP = [relFP_f;relFP_r];
labels = [cellfun(@(x) [x,'_f'], targetRxns,'UniformOutput',0);cellfun(@(x) [x,'_r'], targetRxns,'UniformOutput',0)];

i = 6;
myrxn = EV8_rxns(mismatch_pri(i));
myrxn2 = regexprep(myrxn,'_','');
FP_EV8 = cellfun(@str2num ,([EV8{strcmp(EV8.ID,myrxn2),[9,8,11,14,10,13,12]}]));
FP_XL = rFP(strcmp(labels,myrxn),:);
plot(FP_XL,FP_EV8,'r.')

%==> all the differences are NUMERICAL DIFFERENCES!


