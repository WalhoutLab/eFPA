TsTbl = readtable('input/suppTbls/proteinTSscore.xlsx');
% set NA to 0 (population mean)
tmp = TsTbl{:,5:end};
tmp(isnan(tmp)) = 0;
TsTbl{:,5:end} = tmp;

% % plan A
% expTbl = TsTbl;
% expTbl{:,5:end} = 1.2.^tmp;

% plan B
FcTbl = readtable('input/suppTbls/proteinTissueMedian_log2FC_againts_reference.xlsx');
popMean_RNA = readtable('input/suppTbls/populationMeans.xlsx','Sheet','RNA');
% align the tables
colLabels = FcTbl.Properties.VariableNames(2:end);
rowlabels = intersect(intersect(FcTbl.gene_id,TsTbl.ensembl_id),popMean_RNA.ensembl_id);
FCmat = FcTbl{:,colLabels};
[A B] = ismember(rowlabels,FcTbl.gene_id);
FCmat = FCmat(B(A),:);
TSmat = TsTbl{:,colLabels};
[A B] = ismember(rowlabels,TsTbl.ensembl_id);
TSmat = TSmat(B(A),:);
popMean = readtable('input/suppTbls/populationMeans.xlsx');
[A B] = ismember(rowlabels,popMean.ensembl_id);
popMeanPro = popMean.prt_fitted_mu0(B(A));
[A B] = ismember(rowlabels,popMean_RNA.ensembl_id);
popMeanRNA = popMean_RNA.rna_fitted_mu0(B(A));

% recenter the fc matrix
FCmat_centered = FCmat - popMeanPro;
% set NA to 0 (no FC)
FCmat_centered(isnan(FCmat_centered)) = 0;
% dont weight
FCmat_weighted = 1.*(TSmat~=0) .* FCmat_centered;%(1-2.^(-abs(TSmat))) .* FCmat_centered;
% make the output table
[A B] = ismember(rowlabels, TsTbl.ensembl_id);
expTbl = TsTbl(B(A),[{'ensembl_id','entrez_id','hgnc_name','hgnc_symbol'},colLabels]);
expTbl{:,5:end} = 2.^FCmat_weighted .* (2 .^ popMeanRNA);

expTbl_raw = expTbl;



% ```load expression files```
% load the TS score
TsTbl = readtable('input/suppTbls/proteinTSscore.xlsx');
% set NA to 0 (population mean)
tmp = TsTbl{:,5:end};
tmp(isnan(tmp)) = 0;
TsTbl{:,5:end} = tmp;

% % plan A
% expTbl = TsTbl;
% expTbl{:,5:end} = 1.2.^tmp;

% plan B
% load the median abundance
FcTbl = readtable('input/suppTbls/proteinTissueMedian_log2FC_againts_reference.xlsx');
popMean_RNA = readtable('input/suppTbls/populationMeans.xlsx','Sheet','RNA');
% align the tables
colLabels = FcTbl.Properties.VariableNames(2:end);
rowlabels = intersect(intersect(FcTbl.gene_id,TsTbl.ensembl_id),popMean_RNA.ensembl_id);
FCmat = FcTbl{:,colLabels};
[A B] = ismember(rowlabels,FcTbl.gene_id);
FCmat = FCmat(B(A),:);
TSmat = TsTbl{:,colLabels};
[A B] = ismember(rowlabels,TsTbl.ensembl_id);
TSmat = TSmat(B(A),:);
popMean = readtable('input/suppTbls/populationMeans.xlsx');
[A B] = ismember(rowlabels,popMean.ensembl_id);
popMeanPro = popMean.prt_fitted_mu0(B(A));
[A B] = ismember(rowlabels,popMean_RNA.ensembl_id);
popMeanRNA = popMean_RNA.rna_fitted_mu0(B(A));

% weight the fold change by TS score
FCmat_centered = FCmat - popMeanPro;
% set NA to 0 (no FC)
FCmat_centered(isnan(FCmat_centered)) = 0;
% weight
FCmat_weighted = (2.*normcdf(abs(TSmat))-1) .* FCmat_centered;
% make the output table
[A B] = ismember(rowlabels, TsTbl.ensembl_id);
expTbl = TsTbl(B(A),[{'ensembl_id','entrez_id','hgnc_name','hgnc_symbol'},colLabels]);
expTbl{:,5:end} = 2.^FCmat_weighted .* (2 .^ popMeanRNA);


expTbl_shrunk = expTbl;



%% prepare model
load('input/ihuman_serum.mat');
model.subSystems = [model.subSystems{:}]';
% the default constraints are unlimited, so we dont change anything

% ```special treatments for iHuman1```
% Remove parentathsis in the reaction ID (which will be changed in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% change the name starting with numbers 
model.rxns = regexprep(model.rxns,'(^[0-9])','x$1');


% find the target internal reactions 
exludeSys = {'Transport reactions','Exchange/demand reactions'};
upkrxns = model.rxns(any(model.S(ismember(model.mets,{'majorNutr','sideNutr'}),:),1));
intRxns = model.rxns(~ismember(model.subSystems,exludeSys) & ~ismember(model.rxns,upkrxns));

% create missing field
model = creategrRulesField(model);


targetRxns = intRxns;%{'LEUTA','VALTAm','OIVD1m','OIVD2m','r0655','MCCCrm','MGCHrm','HMGLm'};
conditions = expTbl.Properties.VariableNames(5:end); % conditions to calculate FPA on
