%% load the variables for downstream analysis
% load protein TS score
TsTbl = readtable('input/suppTbls/proteinTSscore.xlsx');
% set NA to 0 (population mean)
% tmp = TsTbl{:,5:end};
% tmp(isnan(tmp)) = 0;
% TsTbl{:,5:end} = tmp;
% load tissue medians
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
% first load the raw FC: dont weight
FCmat_weighted = 1.*(~isnan(TSmat)) .* FCmat_centered;% To ensure it is comparable with the compressed version, we put it to zero when TS score is not available 
% make the output table
[A B] = ismember(rowlabels, TsTbl.ensembl_id);
expTbl = TsTbl(B(A),[{'ensembl_id','entrez_id','hgnc_name','hgnc_symbol'},colLabels]);
expTbl{:,5:end} = 2.^FCmat_weighted .* (2 .^ popMeanRNA);

expTbl_raw = expTbl;

% prepare the compressed tissue medians
% set NA TS score to 0 (population mean)
tmp = TsTbl{:,5:end};
tmp(isnan(tmp)) = 0;
TsTbl{:,5:end} = tmp;
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
%% load the model
load('input/ihuman_serum.mat');
model.subSystems = [model.subSystems{:}]';
% the default constraints are unlimited, so we dont change anything

% ```special treatments for iHuman1```
% Remove parentathsis in the reaction ID (which will be changed in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% change the name starting with numbers 
model.rxns = regexprep(model.rxns,'(^[0-9])','x$1');

% define reaction sets
% the exchange reactions with environment 
excRxns = model.rxns(findExcRxns_XL(model));
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
excRxns = intersect(excRxns,EXinvolvedRxns);

% the transporter (with env) reactions
allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
TSP = [];
for i = 1:length(allCmp_iHumanName)
    myMet_e = {[allCmp_iHumanName{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,allCmp_iHumanName(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    % we define the transporter as the reactions that contain the same
    % metabolite in [e] and another compartment (cellular) in the same
    % reaction
    candidate = intersect(myRxns_non_e,myRxns_e);
    % check if is on diff side of the reaction
    if ~isempty(candidate)
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                TSP = union(TSP,candidate(j));
            end
        end
    end
end
% some special transporter will be missed, we add back 
envTspRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
envTspRxns = intersect(envTspRxns,EXinvolvedRxns);
tspRxns = union(TSP, envTspRxns);

% the internal transporters 
% internal transporters are hard to define, since a reaction can span two
% compartments internally but it is not a real transporter. So, we use the
% subsys annotation as a compromise
intTspRxns = setdiff(model.rxns(ismember(model.subSystems,{'Transport reactions'})),tspRxns);

% internal regular reactions
intRxns = setdiff(model.rxns, [excRxns; tspRxns; intTspRxns]);

% create missing field
model = creategrRulesField(model);

% other variables
targetRxns = intRxns;
conditions = expTbl.Properties.VariableNames(5:end); % conditions to calculate FPA on
