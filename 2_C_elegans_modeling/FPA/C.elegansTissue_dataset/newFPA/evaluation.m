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
load FPA_rxns_oriMERGE.mat
relFP_f_ori = relFP_f;
relFP_r_ori = relFP_r;
load FPA_rxns_wtdDist.mat
relFP_f_wtdDist = relFP_f;
relFP_r_wtdDist = relFP_r;
load FPA_rxns_oriDist_exp_decay_base100.mat
relFP_f_exp_decay100 = relFP_f;
relFP_r_exp_decay100 = relFP_r;
load FPA_rxns_oriDist_exp_decay_base2.mat
relFP_f_exp_decay2 = relFP_f;
relFP_r_exp_decay2 = relFP_r;
load FPA_rxns_wtdDist_exp_decay_base100.mat
relFP_f_wtdDist_exp_decay100 = relFP_f;
relFP_r_wtdDist_exp_decay100 = relFP_r;
load FPA_rxns_wtdDist_exp_decay_base2.mat
relFP_f_wtdDist_exp_decay2 = relFP_f;
relFP_r_wtdDist_exp_decay2 = relFP_r;
%% plot
metainfo = readtable('controlList.xlsx');
qryRxn = 'collg[c]';

info = unique(metainfo.What(strcmp(metainfo.ID,qryRxn)));
if strcmp(qryRxn(end),'r')  
    direction = '_r';
    qryRxn = qryRxn(1:end-1);
else
    direction = '_f';
    qryRxn = regexprep(qryRxn,'f$','');
end

rFP1 = eval(['relFP',direction,'_ori']);
rFP1 = rFP1(strcmp(targetRxns,qryRxn),:);
rFP2 = eval(['relFP',direction,'_wtdDist']);
rFP2 = rFP2(strcmp(targetRxns,qryRxn),:);
rFP3 = eval(['relFP',direction,'_wtdDist_exp_decay2']);
rFP3 = rFP3(strcmp(targetRxns,qryRxn),:);

figure(1)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP1)
title({[qryRxn,direction,' - original']
    [info{1}]})

figure(2)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP2)
title({[qryRxn,direction,' - wtdDist']
    [info{1}]})

figure(3)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP3)
title({[qryRxn,direction,' - wtdDist+expDecay']
    [info{1}]})


%% reactions not in Fig 7
qryRxn = {'TCE0342_f'}


qryRxn = regexprep(qryRxn,'_','');
qryRxn = qryRxn{:};
printRxnFormula(model,[qryRxn(1:end-1),'_X']);
printRxnFormula_XL(model,[qryRxn(1:end-1),'_X']);

if strcmp(qryRxn(end),'r')  
    direction = '_r';
    qryRxn = qryRxn(1:end-1);
else
    direction = '_f';
    qryRxn = regexprep(qryRxn,'f$','');
end

rFP1 = eval(['relFP',direction,'_ori']);
rFP1 = rFP1(strcmp(targetRxns,qryRxn),:);
rFP2 = eval(['relFP',direction,'_wtdDist']);
rFP2 = rFP2(strcmp(targetRxns,qryRxn),:);
rFP3 = eval(['relFP',direction,'_exp_decay2']);
rFP3 = rFP3(strcmp(targetRxns,qryRxn),:);
rFP4 = eval(['relFP',direction,'_exp_decay100']);
rFP4 = rFP4(strcmp(targetRxns,qryRxn),:);
rFP5 = eval(['relFP',direction,'_wtdDist_exp_decay2']);
rFP5 = rFP5(strcmp(targetRxns,qryRxn),:);
rFP6 = eval(['relFP',direction,'_wtdDist_exp_decay100']);
rFP6 = rFP6(strcmp(targetRxns,qryRxn),:);

figure(1)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP1)
title({[qryRxn,direction,' - original']
    })

figure(2)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP2)
title({[qryRxn,direction,' - wtdDist']
    })

figure(3)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP3)
title({[qryRxn,direction,' - expDecay(2)']
    })


figure(4)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP4)
title({[qryRxn,direction,' - expDecay(100)']
    })


figure(5)
c = categorical(regexprep(tissueLabel,'_',' '));
bar(c,rFP5)
title({[qryRxn,direction,' - wtdDist+expDecay(2)']
    })


figure(6)
c = categorical(regexprep(tissueLabel,'_',' '));    
bar(c,rFP6)
title({[qryRxn,direction,' - wtdDist+expDecay(100)']
    })


%% compare the primary and secondary sites
load FPA_interpretation_oriMERGE.mat
labels_ori = labels;
primarySites_ori = primarySites;
SecondSites_ori = SecondSites;
load FPA_interpretation_wtdDist.mat
labels_wtdDist = labels;
primarySites_wtdDist = primarySites;
SecondSites_wtdDist = SecondSites;
load FPA_interpretation_wtdDist_exp_decay.mat
labels_wtdDist_exp_decay = labels;
primarySites_wtdDist_exp_decay = primarySites;
SecondSites_wtdDist_exp_decay = SecondSites;

allrxns = union(union(labels_ori,labels_wtdDist),labels_wtdDist_exp_decay);
primarySitesCmp = repmat({'NA'},length(allrxns), 3);
SecondSitesCmp = repmat({'NA'},length(allrxns), 3);
[A B] = ismember(allrxns,labels_ori);
primarySitesCmp(A,1) = primarySites_ori(B(A));
SecondSitesCmp(A,1) = SecondSites_ori(B(A));
[A B] = ismember(allrxns,labels_wtdDist);
primarySitesCmp(A,2) = primarySites_wtdDist(B(A));
SecondSitesCmp(A,2) = SecondSites_wtdDist(B(A));
[A B] = ismember(allrxns,labels_wtdDist_exp_decay);
primarySitesCmp(A,3) = primarySites_wtdDist_exp_decay(B(A));
SecondSitesCmp(A,3) = SecondSites_wtdDist_exp_decay(B(A));

% clear all NA rxns 
ind1 = all(strcmp(primarySitesCmp,'NA'),2) & all(strcmp(SecondSitesCmp,'NA'),2);
primarySitesCmp(ind1,:) = [];
SecondSitesCmp(ind1,:) = [];
allrxns(ind1) = [];

% save table
mytbl = table(allrxns,primarySitesCmp(:,1),primarySitesCmp(:,2),primarySitesCmp(:,3),SecondSitesCmp(:,1),SecondSitesCmp(:,2),SecondSitesCmp(:,3),...
    'VariableNames',{'rxnID','primarySites_ori','primarySites_wtdDist','primarySites_wtdDist_exp_decay',...
    'SecondSites_ori','SecondSites_wtdDist','SecondSites_wtdDist_exp_decay'});

writetable(mytbl,'interpretationTbl.csv');

% check for some metrics
same_primary = zeros(length(allrxns),2);
same_secondary = zeros(length(allrxns),2);
for i = 1:length(allrxns)
    if strcmp(primarySitesCmp(i,2),primarySitesCmp(i,1))
        same_primary(i,1) = 1;
    end
    if strcmp(primarySitesCmp(i,3),primarySitesCmp(i,1))
        same_primary(i,2) = 1;
    end
    if strcmp(SecondSitesCmp(i,2),SecondSitesCmp(i,1))
        same_secondary(i,1) = 1;
    end
    if strcmp(SecondSitesCmp(i,3),SecondSitesCmp(i,1))
        same_secondary(i,2) = 1;
    end
end
% how many sites are identified
fprintf('original identified %d primary and %d secondary sites\n',sum(~strcmp(primarySitesCmp(:,1),'NA')),sum(~strcmp(SecondSitesCmp(:,1),'NA')));
fprintf('wtdDist identified %d primary and %d secondary sites\n',sum(~strcmp(primarySitesCmp(:,2),'NA')),sum(~strcmp(SecondSitesCmp(:,2),'NA')));
fprintf('wtdDist+expDecay identified %d primary and %d secondary sites\n',sum(~strcmp(primarySitesCmp(:,3),'NA')),sum(~strcmp(SecondSitesCmp(:,3),'NA')));

fprintf('%d original primary sites are recalled in wtdDist\n',sum(~strcmp(primarySitesCmp(:,2),'NA') & same_primary(:,1)));
fprintf('%d original primary sites are recalled in wtdDist+expDecay\n',sum(~strcmp(primarySitesCmp(:,3),'NA') & same_primary(:,2)));
fprintf('%d original secondary sites are recalled in wtdDist\n',sum(~strcmp(SecondSitesCmp(:,2),'NA') & same_secondary(:,1)));
fprintf('%d original secondary sites are recalled in wtdDist+expDecay\n',sum(~strcmp(SecondSitesCmp(:,3),'NA') & same_secondary(:,2)));


fprintf('%d rxns have identical primary and secondary sites in 3 methods\n',sum(all(same_primary,2) & all(same_secondary,2) & ~strcmp(primarySitesCmp(:,1),'NA') & ~strcmp(SecondSitesCmp(:,1),'NA')));
fprintf('%d rxns have identical primary in 3 methods\n',sum(all(same_primary,2)& ~strcmp(primarySitesCmp(:,1),'NA')));
fprintf('%d rxns have identical secondary sites in 3 methods\n',sum(all(same_secondary,2)& ~strcmp(SecondSitesCmp(:,1),'NA')));

% metabolite level
metObj1 = allrxns(cellfun(@(x) ~isempty(regexp(x,'^(UPK|TCE|SNK|DMN)','once')),allrxns));
metObj2 = allrxns(cellfun(@(x) ~isempty(regexp(x,'\[.\]','once')),allrxns)); % targets for metabolite-centric calculation
metInd = ismember(allrxns,[metObj1; metObj2]);

fprintf('original identified %d primary and %d secondary mets\n',sum(metInd & ~strcmp(primarySitesCmp(:,1),'NA')),sum(metInd & ~strcmp(SecondSitesCmp(:,1),'NA')));
fprintf('wtdDist identified %d primary and %d secondary mets\n',sum(metInd & ~strcmp(primarySitesCmp(:,2),'NA')),sum(metInd & ~strcmp(SecondSitesCmp(:,2),'NA')));
fprintf('wtdDist+expDecay identified %d primary and %d secondary mets\n',sum(metInd & ~strcmp(primarySitesCmp(:,3),'NA')),sum(metInd & ~strcmp(SecondSitesCmp(:,3),'NA')));

fprintf('%d original primary mets are recalled in wtdDist\n',sum(metInd & ~strcmp(primarySitesCmp(:,2),'NA') & same_primary(:,1)));
fprintf('%d original primary mets are recalled in wtdDist+expDecay\n',sum(metInd & ~strcmp(primarySitesCmp(:,3),'NA') & same_primary(:,2)));
fprintf('%d original secondary mets are recalled in wtdDist\n',sum(metInd & ~strcmp(SecondSitesCmp(:,2),'NA') & same_secondary(:,1)));
fprintf('%d original secondary mets are recalled in wtdDist+expDecay\n',sum(metInd & ~strcmp(SecondSitesCmp(:,3),'NA') & same_secondary(:,2)));



