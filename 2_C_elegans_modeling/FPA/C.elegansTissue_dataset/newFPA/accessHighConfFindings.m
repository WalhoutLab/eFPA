algorithms = {'oriMERGE','wtdDist','oriDist_exp_decay_base2','oriDist_exp_decay_base100','wtdDist_exp_decay_base2','wtdDist_exp_decay_base100'};
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

%% make the heatmap
% label mismatch??
% first sort by the confidence level for y axis
indexY = [find(strcmp(metainfo.Level,'known'));find(strcmp(metainfo.Level,'coherent'));find(strcmp(metainfo.Level,'unexplored'))];
levels = levels(indexY,:);
yLabels = regexprep(annotationStr(indexY),'_','\\_');
xLabels = regexprep({'oriMERGE','wtdDist','oriDist_exp_decay_base2','oriDist_exp_decay_base100','wtdDist_exp_decay_base2','wtdDist_exp_decay_base100'},'_','\\_');;
% heatmap 
a = heatmap(xLabels,yLabels,levels);
a.Colormap = [0 0 0;
              102/255 204/255 0;
              255/255 255/255 0;];
%% make the dendrigram of pri/sec distance or overall distance



