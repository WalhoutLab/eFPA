algorithms = {'oriMERGE','wtdDist','wtdDist_exp_decay_base2','wtdDist_exp_decay_base100'};
metainfo = readtable('controlList.xlsx');
metainfo.Tissue = regexprep(metainfo.Tissue,'glia','Glia');
metainfo.Tissue = regexprep(metainfo.Tissue,'gonad','Gonad');
metainfo.Tissue = regexprep(metainfo.Tissue,'hypodermis','Hypodermis');
metainfo.Tissue = regexprep(metainfo.Tissue,'intestine','Intestine');
metainfo.Tissue = regexprep(metainfo.Tissue,'muscle','Body_wall_muscle');
metainfo.Tissue = regexprep(metainfo.Tissue,'neurons','Neurons');
metainfo.Tissue = regexprep(metainfo.Tissue,'pharynx','Pharynx');
annotationStr = {};
for i = 1:length(metainfo.ID)
    annotationStr{i,1} = [metainfo.Tissue{i},'-',metainfo.ID{i},' [',metainfo.Where{i},']','[',metainfo.Level{i},']'];
end

levels = zeros(length(metainfo.ID), length(algorithms));
% zero ==> no call; one ==> primary; two ==> secondary 
for i = 1:length(algorithms)
    load(['FPA_interpretation_',algorithms{i},'.mat']);
    labels = regexprep(labels, '_f$','f');
    labels = regexprep(labels, '_r$','r');
    labels = regexprep(labels, '].$',']');
    for j = 1:length(metainfo.ID)
        rxnID = metainfo.ID{j};
        tissue = metainfo.Tissue{j};
        if any(strcmp(rxnID,labels))
            callPri = primarySites{strcmp(rxnID,labels)};
            callSec = SecondSites{strcmp(rxnID,labels)};
        elseif any(strcmp(rxnID,rxns_error))
            fprintf('error seen for rxn %s in algorithm %s\n',rxnID,algorithms{i});
        else
            callPri = 'NA';
            callSec = 'NA';
        end
        if strcmp(tissue,callPri)
            levels(j,i) = 1;
        elseif strcmp(tissue,callSec)
            levels(j,i) = 2;
        end
    end
end
levels_ori = levels;
%% make the heatmap
% first sort by the confidence level for y axis
indexY = [find(strcmp(metainfo.Level,'known'));find(strcmp(metainfo.Level,'coherent'));find(strcmp(metainfo.Level,'unexplored'))];
levels = levels_ori(indexY,:);
yLabels = regexprep(annotationStr(indexY),'_','\\_');
xLabels = regexprep({'original','original + weighted distance','improved FPA (base = 2)','improved FPA (base = 100)'},'_','\\_');
% heatmap
levels4cluster = levels;
levels4cluster(levels4cluster==1) = 3;
Y = pdist(levels4cluster);
Z = linkage(Y);
figure;
[~,~,labelClustered] = dendrogram(Z,1000);
figure('units','inch','position',[0,0,8,16])
a = heatmap(xLabels,yLabels(labelClustered),levels(labelClustered,:),'ColorbarVisible',0,'CellLabelColor','none');
a.Colormap = [0 0 0;
              102/255 204/255 0;
              255/255 255/255 0;];
a.FontSize = 10;
hAx=a.NodeChildren(3);
hAx.XAxis.FontWeight='bold';
hAx.XAxis.FontSize=12;
S = hgexport('readstyle','default_sci');
style.Format = 'tiff';
style.ApplyStyle = '1';
hgexport(gcf,'test',S,'applystyle',true);
saveas(gca,'figures/benchmarkNewFPA_heatmap.tiff');
%% make the dendrigram of pri/sec distance or overall distance
tblEV7 = readtable('Table_EV7.xlsx');
tblEV8 = readtable('Table_EV8.xlsx');
allCalls = [tblEV7(:,[1 5 6]); tblEV8(:,[1 6 7])];
starred = allCalls(contains(allCalls.Primary,'*') | contains(allCalls.Secondary,'*'),:);
refCalls = starred(:,1:2);
refCalls.Properties.VariableNames = {'ID','Secondary'};
refCalls = [refCalls; starred(:,[1 3])];
refCalls.Properties.VariableNames = {'ID','Tissue'};
refCalls(~contains(refCalls.Tissue,'*'),:) = [];
refCalls.Tissue = regexprep(refCalls.Tissue,'\*$','');
refCalls.Tissue = regexprep(refCalls.Tissue,'gli','Glia');
refCalls.Tissue = regexprep(refCalls.Tissue,'gon','Gonad');
refCalls.Tissue = regexprep(refCalls.Tissue,'hyp','Hypodermis');
refCalls.Tissue = regexprep(refCalls.Tissue,'int','Intestine');
refCalls.Tissue = regexprep(refCalls.Tissue,'mus','Body_wall_muscle');
refCalls.Tissue = regexprep(refCalls.Tissue,'neu','Neurons');
refCalls.Tissue = regexprep(refCalls.Tissue,'pha','Pharynx');
%% only the distance for primary/secondary sites in original FPA
levels = zeros(length(metainfo.ID), length(algorithms));
% zero ==> no call; one ==> primary; two ==> secondary 
for i = 1:length(algorithms)
    load(['FPA_interpretation_',algorithms{i},'.mat']);
    labels = regexprep(labels, '_f$','f');
    labels = regexprep(labels, '_r$','r');
    labels = regexprep(labels, '].$',']');
    for j = 1:length(refCalls.ID)
        rxnID = refCalls.ID{j};
        tissue = refCalls.Tissue{j};
        if any(strcmp(rxnID,labels))
            callPri = primarySites{strcmp(rxnID,labels)};
            callSec = SecondSites{strcmp(rxnID,labels)};
        elseif any(strcmp(rxnID,rxns_error))
            fprintf('error seen for rxn %s in algorithm %s\n',rxnID,algorithms{i});
        else
            callPri = 'NA';
            callSec = 'NA';
        end
        if strcmp(tissue,callPri) || strcmp(tissue,callSec)
            levels(j,i) = 1;
%         elseif strcmp(tissue,callSec)
%             levels(j,i) = 2;
        end
    end
end
levels_ori = levels;

Y = pdist(levels','jaccard');% note: 0 is not counted in Jaccard, so that it wouldnt be biased by the unbalanced large number of no calls 
Z = linkage(Y);
figure;
[h,~,labelClustered] = dendrogram(Z,10000,'Labels' ,{'original','original + weighted distance','improved FPA (base = 2)','improved FPA (base = 100)'});
xtickangle(30)
h(1).Color = 'k';
h(2).Color = 'k';
h(3).Color = 'k';
ylabel('Jaccard distance');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/benchmarkNewFPA_dendr_by_starSites_in_MERGE.tiff']);
%% the distance for all rxn's primary/secondary/none site calls 
% we only distinguish pri/sec call from no calls
allTissueCalls = zeros(size(allCalls,1)*7,1);
tissues = {'gli','gon','hyp','int','mus','neu','pha'};
priCalls = regexprep(allCalls.Primary,'\*$','');
SecCalls = regexprep(allCalls.Secondary,'\*$','');
for i = 1:length(tissues)
    isSite = find(strcmp(priCalls,tissues{i}) | strcmp(SecCalls,tissues{i})); 
    allTissueCalls((i-1) * size(allCalls,1) + isSite) = 1;
end
allTissueCalls_ID = repmat(allCalls.ID,7,1);
allTissueCalls_tissue = [repmat({'Glia'},size(allCalls,1),1);
                         repmat({'Gonad'},size(allCalls,1),1);
                         repmat({'Hypodermis'},size(allCalls,1),1);
                         repmat({'Intestine'},size(allCalls,1),1);
                         repmat({'Body_wall_muscle'},size(allCalls,1),1);
                         repmat({'Neurons'},size(allCalls,1),1);
                         repmat({'Pharynx'},size(allCalls,1),1)];

levels = zeros(length(allTissueCalls), length(algorithms));
% zero ==> no call; one ==> primary; two ==> secondary 
for i = 1:length(algorithms)
    load(['FPA_interpretation_',algorithms{i},'.mat']);
    labels = regexprep(labels, '_f$','f');
    labels = regexprep(labels, '_r$','r');
    labels = regexprep(labels, '].$',']');
    for j = 1:length(allTissueCalls_ID)
        rxnID = allTissueCalls_ID{j};
        tissue = allTissueCalls_tissue{j};
        if any(strcmp(rxnID,labels))
            callPri = primarySites{strcmp(rxnID,labels)};
            callSec = SecondSites{strcmp(rxnID,labels)};
        elseif any(strcmp(rxnID,rxns_error))
            fprintf('error seen for rxn %s in algorithm %s\n',rxnID,algorithms{i});
        else
            callPri = 'NA';
            callSec = 'NA';
        end
        if strcmp(tissue,callPri) || strcmp(tissue,callSec)
            levels(j,i) = 1;
%         elseif strcmp(tissue,callSec)
%             levels(j,i) = 2;
        end
    end
end
levels_ori = levels;
% check the distance between replicate of MERGE and Safak's MERGE
pdist([allTissueCalls,levels(:,1)]','jaccard')
Y = pdist(levels','jaccard');% note: 0 is not counted in Jaccard, so that it wouldnt be biased by the unbalanced large number of no calls 
Z = linkage(Y);
figure;
[h,~,labelClustered] = dendrogram(Z,10000,'Labels' ,{'original','original + weighted distance','improved FPA (base = 2)','improved FPA (base = 100)'});
xtickangle(30)
h(1).Color = 'k';
h(2).Color = 'k';
h(3).Color = 'k';
ylabel('Jaccard distance');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/benchmarkNewFPA_dendr_by_allSites.tiff']);
%% overall distance for rFP values 
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
ID2 = allCalls.ID;
[A B] = ismember(ID2,targetMets);
ID2(A) = cellfun(@(x) [x,'f'], ID2(A),'UniformOutput',0);
[A B] = ismember(ID2, labels);
levels = zeros(length(allCalls.ID(A))*7, length(algorithms));

for i = 1:length(algorithms)
    load(['FPA_rxns_',algorithms{i},'.mat']);
    rFP = [relFP_f;relFP_r];
    for j = 1:7
        levels(((j-1)*length(allCalls.ID(A)) + 1) : j*length(allCalls.ID(A)),i) = rFP(B(A),j);
    end
end

% calculate cosine distance
levels(all(isnan(levels),2),:) = [];
Y = pdist(levels','cosine');% note: 0 is not counted in Jaccard, so that it wouldnt be biased by the unbalanced large number of no calls 
Z = linkage(Y);
figure;
[h,~,labelClustered] = dendrogram(Z,10000,'Labels' ,{'original','original + weighted distance','improved FPA (base = 2)','improved FPA (base = 100)'});
xtickangle(30)
h(1).Color = 'k';
h(2).Color = 'k';
h(3).Color = 'k';
ylabel('Cosine distance');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.2, 4.3];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/benchmarkNewFPA_dendr_by_rFP_cosine.tiff']);
 
    
    
    
    
