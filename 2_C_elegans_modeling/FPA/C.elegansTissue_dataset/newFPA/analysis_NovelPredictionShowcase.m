%% load the rFP data
algorithms = {'oriMERGE','wtdDist','wtdDist_exp_decay_base2','wtdDist_exp_decay_base100'};
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

labels = [cellfun(@(x) [x,'f'], targetRxns,'UniformOutput',0);cellfun(@(x) [x,'r'], targetRxns,'UniformOutput',0)];
rFPs = struct();
for i = 1:length(algorithms)
    load(['FPA_rxns_',algorithms{i},'.mat']);
    rFPs.(algorithms{i}) = [relFP_f;relFP_r];
end

%% plot case1: production of pail3p
% https://www.frontiersin.org/articles/10.3389/fnmol.2019.00208/full
tissueLabel2 = regexprep(tissueLabel,'_',' ');
order0 = [1 2 6 7 4 5 3];
qryRxn = 'DMN0054f';
figure;
hold on
myLevel = rFPs.oriMERGE(strcmp(labels,qryRxn),:);
X = categorical(tissueLabel2(order0));
X = reordercats(X,tissueLabel2(order0));
Y1 = myLevel(order0);
% other algorithms
myLevel = rFPs.wtdDist(strcmp(labels,qryRxn),:);
Y2 = myLevel(order0);
myLevel = rFPs.wtdDist_exp_decay_base100(strcmp(labels,qryRxn),:);
Y3 = myLevel(order0);
myLevel = rFPs.wtdDist_exp_decay_base2(strcmp(labels,qryRxn),:);
Y4 = myLevel(order0);
% bars
b = bar(X,[Y1',Y2',Y3',Y4']);
b(1).FaceColor = '#0072BD';
b(2).FaceColor = '#EDB120';
b(3).FaceColor = '#D95319';
b(4).FaceColor = '#77AC30';
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [13, 8];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
legend({'original','original + weighted distance','improved FPA (base = 2)','improved FPA (base = 100)'},'FontSize',12)
ylabel('Relative flux potential',  'FontSize',15)
plt.Interpreter = 'none';
plt.LegendLoc = 'northeast';
plt.export('figures/NovelPred_pail3p.tiff');

%% display the expression of genes in the pathway 
close all
geneset = {'ttx-7','gpap-1','cdgs-1','vps-34','mtm-3', 'pisy-1','inos-1'};% 
for i = 1:length(geneset)
    qryGene = geneset{i};
    order0 = [1 2 6 7 4 5 3];
    figure;
    hold on
    myLevel = expressionTbl{strcmp(expressionTbl.Properties.RowNames,qryGene),:};
    Y1 = myLevel(order0);
    % bars
    b = bar(X,Y1);
    title(qryGene);
    b(1).FaceColor = [0.5 0.5 0.5];
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [3.6, 3];
    plt.LineWidth = 2;
    plt.FontSize = 12;
    plt.FontName = 'Arial';
    ylabel('TPM',  'FontSize',12)
    plt.Interpreter = 'none';
    plt.export(['figures/geneExp/',qryGene,'.tiff']);
end


