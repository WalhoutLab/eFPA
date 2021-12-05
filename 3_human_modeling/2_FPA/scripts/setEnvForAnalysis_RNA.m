TsTbl = readtable('input/suppTbls/RNATSscore.xlsx');
% set NA to 0 (population mean)
tmp = TsTbl{:,5:end};
tmp(isnan(tmp)) = 0;
TsTbl{:,5:end} = tmp;

% % plan A
% expTbl = TsTbl;
% expTbl{:,5:end} = 1.2.^tmp;

% plan B
% load the median abundance
logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
popMean = readtable('input/suppTbls/populationMeans.xlsx','Sheet','RNA');
% align the tables
colLabels = logTPMTbl.Properties.VariableNames(2:end);
rowlabels = intersect(intersect(logTPMTbl.gene_id,TsTbl.ensembl_id),popMean.ensembl_id);
logTPMmat = logTPMTbl{:,colLabels};
[A B] = ismember(rowlabels,logTPMTbl.gene_id);
logTPMmat = logTPMmat(B(A),:);
TSmat = TsTbl{:,colLabels};
[A B] = ismember(rowlabels,TsTbl.ensembl_id);
TSmat = TSmat(B(A),:);
[A B] = ismember(rowlabels,popMean.ensembl_id);
popMeanRNA = popMean.rna_fitted_mu0(B(A));

% weight the fold change by TS score
% FCmat = logTPMmat; %- popMeanRNA;
% set NA to 0 (no FC)
% FCmat(isnan(FCmat)) = 0;
% weight
% FCmat_weighted = (1-2.^(-abs(TSmat))) .* FCmat;
% make the output table
[A B] = ismember(rowlabels, TsTbl.ensembl_id);
expTbl = TsTbl(B(A),[{'ensembl_id','entrez_id','hgnc_name','hgnc_symbol'},colLabels]);
% make the hypothetical count based on shrunk FC
expTbl{:,5:end} = 2.^(logTPMmat);

expTbl_raw = expTbl;



% ```load expression files```
% load the TS score
TsTbl = readtable('input/suppTbls/RNATSscore.xlsx');
% set NA to 0 (population mean)
tmp = TsTbl{:,5:end};
tmp(isnan(tmp)) = 0;
TsTbl{:,5:end} = tmp;

% % plan A
% expTbl = TsTbl;
% expTbl{:,5:end} = 1.2.^tmp;

% plan B
% load the median abundance
logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
popMean = readtable('input/suppTbls/populationMeans.xlsx','Sheet','RNA');
% align the tables
colLabels = logTPMTbl.Properties.VariableNames(2:end);
rowlabels = intersect(intersect(logTPMTbl.gene_id,TsTbl.ensembl_id),popMean.ensembl_id);
logTPMmat = logTPMTbl{:,colLabels};
[A B] = ismember(rowlabels,logTPMTbl.gene_id);
logTPMmat = logTPMmat(B(A),:);
TSmat = TsTbl{:,colLabels};
[A B] = ismember(rowlabels,TsTbl.ensembl_id);
TSmat = TSmat(B(A),:);
[A B] = ismember(rowlabels,popMean.ensembl_id);
popMeanRNA = popMean.rna_fitted_mu0(B(A));

% weight the fold change by TS score
FCmat = logTPMmat - popMeanRNA;
% set NA to 0 (no FC)
FCmat(isnan(FCmat)) = 0;
% weight
FCmat_weighted = (2.*normcdf(abs(TSmat))-1) .* FCmat;
% make the output table
[A B] = ismember(rowlabels, TsTbl.ensembl_id);
expTbl = TsTbl(B(A),[{'ensembl_id','entrez_id','hgnc_name','hgnc_symbol'},colLabels]);
% make the hypothetical count based on shrunk FC
expTbl{:,5:end} = 2.^(popMeanRNA+FCmat_weighted);

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
intRxns = model.rxns(~ismember(model.subSystems,exludeSys));

% create missing field
model = creategrRulesField(model);


targetRxns = intRxns;%{'LEUTA','VALTAm','OIVD1m','OIVD2m','r0655','MCCCrm','MGCHrm','HMGLm'};
conditions = expTbl.Properties.VariableNames(5:end); % conditions to calculate FPA on
