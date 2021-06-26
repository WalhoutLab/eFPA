%%
for i = 1:32
    sampleName = conditions{i};
    load(['output/humanModel/',lib,'/OFD/',sampleName,'.mat']);
    fprintf('lactate production of %s is %f \n',sampleName,myCSM.OFD(strcmp(model.rxns,'HMR_9135')));
    fprintf('lactate production ratio of %s is %f \n',sampleName,myCSM.OFD(strcmp(model.rxns,'HMR_9135'))/myCSM.OFD(strcmp(model.rxns,'EX_majorNutr')));
end


%%
OFD_mat = [];
for i = 1:32
    sampleName = conditions{i};
    load(['output/humanModel/',lib,'/OFD/',sampleName,'.mat']);
    OFD_mat = [OFD_mat,myCSM.OFD];
end

%%
rowLabels = model.rxns;
rowLabels(var(OFD_mat,[],2)==0,:) = [];
OFD_mat(var(OFD_mat,[],2)==0,:) = [];
OFD_mat = OFD_mat ./ max(abs(OFD_mat),[],2);


cgo=clustergram(OFD_mat,'ColumnLabels',conditions,'RowLabels',rowLabels);
c=get(cgo,'ColorMap');
n = 100;
tmp = [ones(n,1), linspace(1,0,n)',linspace(1,0,n)'];
tmp = tmp(2:end,:);
cpr=[linspace(0,1,n)',linspace(0,1,n)',ones(n,1);...
    tmp];
set(cgo,'ColorMap',cpr);
set(cgo,'Symmetric',true);

set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 12)
%%
ROI = mygroup.RowNodeNames;
ROI = regexprep(ROI,'_.$','');
[ROI_annotated, ROI_enrichment] = annotateRxnSet(ROI,model);