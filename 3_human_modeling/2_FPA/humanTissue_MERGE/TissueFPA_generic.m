function TissueFPA_generic(distanceMat,FPAtype,expressionDataType, distanceOrder, networkType, targetType)
n =  str2num(distanceOrder);
%% This is to reproduce the C. elegans tissue FPA results in supp table S7 and S8
%% set up the env variables
addpath ~/cobratoolbox/
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
initCobraToolbox(false)
mkdir('output');
%% prepare model
load('input/ihuman_serum.mat','model');
model.subSystems = [model.subSystems{:}]';
% the default constraints are unlimited, so we dont change anything

% ```special treatments for iHuman1```
% Remove parentathsis in the reaction ID (which will be changed in distance calculation)
model.rxns = regexprep(model.rxns,'\(|\)|\[|\]|-','_');
% change the name starting with numbers 
model.rxns = regexprep(model.rxns,'(^[0-9])','x$1');

% correct GPR
model = changeGeneAssociation(model, 'HMR_4137',...
    'ENSG00000091140 and ENSG00000110435 and (ENSG00000131828 or ENSG00000163114) and ENSG00000150768 and ENSG00000168291');


% find the target internal reactions 
exludeSys = {'Transport reactions','Exchange/demand reactions'};
upkrxns = model.rxns(any(model.S(ismember(model.mets,{'majorNutr','sideNutr'}),:),1));
intRxns = model.rxns(~ismember(model.subSystems,exludeSys) & ~ismember(model.rxns,upkrxns));


% create missing field
model = creategrRulesField(model);
%% other inputs
% ```load the distance matrix```
if strcmp(distanceMat,'original')
    distance_raw = readtable('./input/distanceMatrix.txt','FileType','text','ReadRowNames',true); % we load from the output of the distance calculator. For usage of distance calculator, please refer to the MetabolicDistance folder
elseif strcmp(distanceMat,'weighted')
    distance_raw = readtable('./input/distanceMatrix_weighted.txt','FileType','text','ReadRowNames',true); % we load from the output of the distance calculator. For usage of distance calculator, please refer to the MetabolicDistance folder
    % there are some NaN. seems related to transport rxns; unknown why, but
    % manually checking a few indicates those reactions are very distal
    % (eg, one is drug metabolism one is ER lipid transporter). so, as a
    % workaround, we set all these distances to inf
    Dmat = table2array(distance_raw);
    Dmat(isnan(Dmat)) = inf;
    distance_raw{:,:} = Dmat;
    clear Dmat
else
    error('unknown distMat type!');
end
labels = distance_raw.Properties.VariableNames;
labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
distMat_raw = table2array(distance_raw);
distMat_min = zeros(size(distance_raw,1),size(distance_raw,2));
for i = 1:size(distMat_min,1)
    for j = 1:size(distMat_min,2)
        distMat_min(i,j) = min([distMat_raw(i,j),distMat_raw(j,i)]);
    end
end
distMat = distMat_min;



% note that the distance matrix dont have distances for uptake reactions as 
% the distance was calculated by the original human1 model. but there is no 
% need to add distances for uptake reactions since their penalty is zero

% ```set the special penalties (if desired)```
extRxns = model.rxns(findExcRxns_XL(model));
manualPenalty = [extRxns];
manualPenalty(:,2) = [mat2cell(zeros(length(extRxns),1),ones(length(extRxns),1))];
manualDist = {};
%% ```load expression files```
if strcmp(expressionDataType,'RNA_TS_common')
    % only look at common genes between RNA and protein dataset 
    % the lowly expressed gene in RNA is also filtered due to lack of pop
    % mean estimation
     
    % load the TS score
    TsTbl = readtable('input/suppTbls/RNATSscore.xlsx');
    TsTbl_pro = readtable('input/suppTbls/proteinTSscore.xlsx');
    % set NA to 0 (population mean)
    tmp = TsTbl{:,5:end};
    tmp(isnan(tmp)) = 0;
    TsTbl{:,5:end} = tmp;
    % load the median abundance
    logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
    popMean = readtable('input/suppTbls/populationMeans.xlsx','Sheet','RNA');
    % align the tables
    colLabels = logTPMTbl.Properties.VariableNames(2:end);
    rowlabels = intersect(intersect(intersect(logTPMTbl.gene_id,TsTbl.ensembl_id),popMean.ensembl_id),TsTbl_pro.ensembl_id);
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
    
elseif strcmp(expressionDataType,'RNA_TS_all') 
    % load the TS score
    TsTbl = readtable('input/suppTbls/RNATSscore.xlsx');
    TsTbl_homemade = readtable('input/suppTbls/RNATS_score_homemade.xlsx');
    % combine the tables
    % first align the colnames 
    TsTbl_homemade.Properties.VariableNames(strcmp(TsTbl_homemade.Properties.VariableNames,'GEjunction')) = {'GEJunction'};
    TsTbl_homemade = TsTbl_homemade(:,[TsTbl_homemade.Properties.VariableNames(1),TsTbl.Properties.VariableNames(5:end)]);
    % then merge
    TsTbl_homemade = TsTbl_homemade(~ismember(TsTbl_homemade.Row,TsTbl.ensembl_id),:);
    TsTbl_homemade.Properties.VariableNames(1) = {'ensembl_id'};
    TsTbl_homemade.entrez_id = repmat(NaN,size(TsTbl_homemade,1),1);
    TsTbl_homemade.hgnc_name = repmat({'Not Determined'},size(TsTbl_homemade,1),1);
    TsTbl_homemade.hgnc_symbol = repmat({'Not Determined'},size(TsTbl_homemade,1),1);
    TsTbl_homemade = TsTbl_homemade(:,TsTbl.Properties.VariableNames);
    TsTbl = [TsTbl;TsTbl_homemade];
    
    % set NA to 0 (population mean)
    tmp = TsTbl{:,5:end};
    tmp(isnan(tmp)) = 0;
    TsTbl{:,5:end} = tmp;
    % load the median abundance
    logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
    popMean = readtable('input/suppTbls/populationMeans.xlsx','Sheet','RNA');
    popMean_homemade = readtable('input/suppTbls/RNATS_fitting_homemade.xlsx');
    popMean = popMean(:,[{'ensembl_id'},{'rna_fitted_mu0'},{'rna_fitted_sd0'}]);
    popMean_homemade = popMean_homemade(:,2:4);
    popMean_homemade.Properties.VariableNames = popMean.Properties.VariableNames;
    popMean = [popMean;popMean_homemade(~ismember(popMean_homemade.ensembl_id, popMean.ensembl_id),:)];
    
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
    
elseif strcmp(expressionDataType,'RNA_raw_common')
    % load the TS score
    TsTbl = readtable('input/suppTbls/RNATSscore.xlsx');
    TsTbl_pro = readtable('input/suppTbls/proteinTSscore.xlsx');
    % set NA to 0 (population mean)
    tmp = TsTbl{:,5:end};
    tmp(isnan(tmp)) = 0;
    TsTbl{:,5:end} = tmp;
    % load the median abundance
    logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
    popMean = readtable('input/suppTbls/populationMeans.xlsx','Sheet','RNA');
    % align the tables
    colLabels = logTPMTbl.Properties.VariableNames(2:end);
    rowlabels = intersect(intersect(intersect(logTPMTbl.gene_id,TsTbl.ensembl_id),popMean.ensembl_id),TsTbl_pro.ensembl_id);
    logTPMmat = logTPMTbl{:,colLabels};
    [A B] = ismember(rowlabels,logTPMTbl.gene_id);
    logTPMmat = logTPMmat(B(A),:);
%     TSmat = TsTbl{:,colLabels};
%     [A B] = ismember(rowlabels,TsTbl.ensembl_id);
%     TSmat = TSmat(B(A),:);
%     [A B] = ismember(rowlabels,popMean.ensembl_id);
    % popMeanRNA = popMean.rna_fitted_mu0(B(A));
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
    
elseif strcmp(expressionDataType,'RNA_raw_all')
    % load the median abundance
    logTPMTbl = readtable('input/suppTbls/RNATissueMedian_log2TPM.xlsx');
    % align the tables
    colLabels = logTPMTbl.Properties.VariableNames(2:end);
    rowlabels = logTPMTbl.gene_id;
    expTbl = table(rowlabels,rowlabels,rowlabels,rowlabels);
    expTbl{:,5:length(colLabels)+4} = 2.^(logTPMTbl{:,2:end});
    expTbl.Properties.VariableNames = [{'ensembl_id','ensembl_id1','ensembl_id2','ensembl_id3'},colLabels];
    
elseif strcmp(expressionDataType,'protein_TS_common')
    % load the TS score
    TsTbl = readtable('input/suppTbls/proteinTSscore.xlsx');
    % set NA to 0 (population mean)
    tmp = TsTbl{:,5:end};
    tmp(isnan(tmp)) = 0;
    TsTbl{:,5:end} = tmp;
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

elseif strcmp(expressionDataType,'protein_raw_common')
    % load the TS score
    TsTbl = readtable('input/suppTbls/proteinTSscore.xlsx');
    % set NA to 0 (population mean)
    tmp = TsTbl{:,5:end};
    tmp(isnan(tmp)) = 0;
    TsTbl{:,5:end} = tmp;
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

elseif strcmp(expressionDataType,'protein_TS_all')
    % load the TS score
    TsTbl = readtable('input/suppTbls/proteinTSscore.xlsx');
    % set NA to 0 (population mean)
    tmp = TsTbl{:,5:end};
    tmp(isnan(tmp)) = 0;
    TsTbl{:,5:end} = tmp;
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

    % assume the expression of undetected RNA is 1 TPM
    undet = setdiff(TsTbl.ensembl_id, popMean_RNA.ensembl_id);
    rowlabels = [rowlabels;undet];
    [A B] = ismember(undet,FcTbl.gene_id);
    FCmat = [FCmat;FcTbl{B(A),colLabels}];
    [A B] = ismember(undet,TsTbl.ensembl_id);
    tmp = TsTbl{:,colLabels};
    TSmat = [TSmat;tmp(B(A),:)];
    [A B] = ismember(undet,popMean.ensembl_id);
    popMeanPro = [popMeanPro;popMean.prt_fitted_mu0(B(A))];
    popMeanRNA = [popMeanRNA;zeros(length(undet),1)];% assume the expression of undetected RNA is 1 TPM
    
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

elseif strcmp(expressionDataType,'protein_raw_all')
    % load the TS score
    TsTbl = readtable('input/suppTbls/proteinTSscore.xlsx');
    % set NA to 0 (population mean)
    tmp = TsTbl{:,5:end};
    tmp(isnan(tmp)) = 0;
    TsTbl{:,5:end} = tmp;
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

    % assume the expression of undetected RNA is 1 TPM
    undet = setdiff(TsTbl.ensembl_id, popMean_RNA.ensembl_id);
    rowlabels = [rowlabels;undet];
    [A B] = ismember(undet,FcTbl.gene_id);
    FCmat = [FCmat;FcTbl{B(A),colLabels}];
    [A B] = ismember(undet,TsTbl.ensembl_id);
    tmp = TsTbl{:,colLabels};
    TSmat = [TSmat;tmp(B(A),:)];
    [A B] = ismember(undet,popMean.ensembl_id);
    popMeanPro = [popMeanPro;popMean.prt_fitted_mu0(B(A))];
    popMeanRNA = [popMeanRNA;zeros(length(undet),1)];% assume the expression of undetected RNA is 1 TPM

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
else
    error('unknown expression data type!');
end

 

% Make the master_expression for human tissues
master_expression = {};
conditions = expTbl.Properties.VariableNames(5:end); % conditions to calculate FPA on
targetGenes = model.genes(ismember(model.genes, expTbl.ensembl_id));
% Fill in the "master_expression" cell array
for i = 1:length(conditions) 
    expression = struct();
    expression.genes = targetGenes; 
    [A B] = ismember(targetGenes,expTbl.ensembl_id);
    expression.value = expTbl.(conditions{i})(B(A));
    master_expression{i} = expression;
end
%% ```setup the block list (context-specific network)```
if strcmp(networkType,'tissue')
    expDataID = strsplit(expressionDataType,'_');
    expDataID = expDataID{1};
    % First, let's merge the level tables for each condition
    for i = 1:length(conditions)
        load(['./../../1_iMAT++/output/humanModel/',expDataID,'/FAA/',conditions{i},'_levels.mat']); % OR use FAA network
        levelTbl_f(:,i) = levels_f;
        levelTbl_r(:,i) = levels_r;
    end
    % Because the FPA is done using the irreversible model, the rxnID is
    % different from the original. We provided a function to get the new rxnIDs
    % to block according to the level tables.
    blockList = getBlockList(model,levelTbl_f,levelTbl_r);
elseif strcmp(networkType,'naive')
    blockList = {};
else
    error('unknown network type!');
end
%% prepare paremeters of running FPA
maxDist = max(distMat(~isinf(distMat))); % maximum distance
%% PART I: FPA for normal reactions (reaction-centered FPA)
if strcmp(targetType,'rxn')
    load('lowFluxRxns.mat');
    load('FVA.mat');
    FVA_label = model.rxns;
    targetRxns = intRxns;%

    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    % FPA with normal alpha
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP] = FPA_clusterWrapper(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList);
    elseif strcmp(FPAtype,'new')
        [FP] = FPA2_oriMERGE_base2_clusterWrapper(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList);
    end

    % FPA for low flux rxns
    changeCobraSolverParams('LP','optTol', 10e-9);
    changeCobraSolverParams('LP','feasTol', 10e-9);
    smallTarget = intersect(lowFluxRxns,targetRxns);
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP_small] = FPA_clusterWrapper(model,smallTarget,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList,...
                                                    {},1,{},0.005);
    elseif strcmp(FPAtype,'new')
        [FP_small] = FPA2_oriMERGE_base2_clusterWrapper(model,smallTarget,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList,...
                                                    {},1,{},0.005);
    end
elseif strcmp(targetType,'transporter')
    load('lowFluxRxns.mat');
    load('FVA.mat');
    FVA_label = model.rxns;
    targetRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
    metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
    metComp = [metComp{:}]';
    EXmets = strcmp(metComp,'[Extracellular]');
    EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
    targetRxns = intersect(targetRxns,EXinvolvedRxns);
    targetRxns = [targetRxns;{'MI1Pt'}];
    % the in-organic co-transportable metabilites are:
    helperMet = {'CO2';'AMP';'NADP+';'NADPH';'PPi';'O2';'NADH';'NAD+';
                'Pi';'ADP';'CoA';'ATP';'H2O';'H+';'GTP';'GDP';
                'Electron Transfer Flavoprotein Reduced';'Electron Transfer Flavoprotein Oxidized';
                'L-carnitine';'FAD';'FADH2';'Na+';'hydroxide';'chloride';'iodide'};
            
    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    if strcmp(FPAtype,'original')
        [FP] = MET_FPA_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    elseif strcmp(FPAtype,'new')
        [FP] = MET_FPA2_oriMERGE_base2_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    end

    % FPA for low flux rxns
    changeCobraSolverParams('LP','optTol', 10e-9);
    changeCobraSolverParams('LP','feasTol', 10e-9);
    smallTarget = intersect(lowFluxRxns,targetRxns);
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP_small] = MET_FPA_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    elseif strcmp(FPAtype,'new')
        [FP_small] = MET_FPA2_oriMERGE_base2_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    end
    
elseif strcmp(targetType,'demand')
    load('allCmp_iHumanName.mat');
    load('FVA_cmp.mat');
    load('lowFluxCmp.mat');
    allCmp_iHumanName = unique(allCmp_iHumanName);
    targetRxns = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);
    FVA_label = targetRxns;
    lowFluxRxns = cellfun(@(x) ['NewMet_',x],lowFluxCmps,'UniformOutput',false);
    % the in-organic co-transportable metabilites are:
    helperMet = {'CO2';'AMP';'NADP+';'NADPH';'PPi';'O2';'NADH';'NAD+';
                'Pi';'ADP';'CoA';'ATP';'H2O';'H+';'GTP';'GDP';
                'Electron Transfer Flavoprotein Reduced';'Electron Transfer Flavoprotein Oxidized';
                'L-carnitine';'FAD';'FADH2';'Na+';'hydroxide';'chloride';'iodide'};
            
    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    if strcmp(FPAtype,'original')
        [FP] = MET_FPA_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    elseif strcmp(FPAtype,'new')
        [FP] = MET_FPA2_oriMERGE_base2_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    end

    % FPA for low flux rxns
    changeCobraSolverParams('LP','optTol', 10e-9);
    changeCobraSolverParams('LP','feasTol', 10e-9);
    smallTarget = intersect(lowFluxRxns,targetRxns);
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP_small] = MET_FPA_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    elseif strcmp(FPAtype,'new')
        [FP_small] = MET_FPA2_oriMERGE_base2_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    end
elseif strcmp(targetType,'allMetDemand')
    load('FVA_cmp_all.mat');
    load('lowFluxCmp_all.mat');
    allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
    allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
    targetRxns = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);
    FVA_label = targetRxns;
    lowFluxRxns = cellfun(@(x) ['NewMet_',x],lowFluxCmps,'UniformOutput',false);
    % the in-organic co-transportable metabilites are:
    helperMet = {'CO2';'AMP';'NADP+';'NADPH';'PPi';'O2';'NADH';'NAD+';
                'Pi';'ADP';'CoA';'ATP';'H2O';'H+';'GTP';'GDP';
                'Electron Transfer Flavoprotein Reduced';'Electron Transfer Flavoprotein Oxidized';
                'L-carnitine';'FAD';'FADH2';'Na+';'hydroxide';'chloride';'iodide'};
            
    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    if strcmp(FPAtype,'original')
        [FP] = MET_FPA_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    elseif strcmp(FPAtype,'new')
        [FP] = MET_FPA2_oriMERGE_base2_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    end

    % FPA for low flux rxns ==> to be done (need to do FVA first) 
    smallTarget = {};
    changeCobraSolverParams('LP','optTol', 10e-9);
    changeCobraSolverParams('LP','feasTol', 10e-9);
    smallTarget = intersect(lowFluxRxns,targetRxns);
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP_small] = MET_FPA_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.0003);
    elseif strcmp(FPAtype,'new')
        [FP_small] = MET_FPA2_oriMERGE_base2_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.0003);
    end
elseif strcmp(targetType,'HMDB')
    output_collections = struct();
    % first finish the rxn-centric task
    load('lowFluxRxns.mat');
    load('FVA.mat');
    FVA_label = model.rxns;
    load('HMDB_objectives.mat');
    targetRxns = targetRxns_INT;
    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    % FPA with normal alpha
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP] = FPA_clusterWrapper(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList);
    elseif strcmp(FPAtype,'new')
        [FP] = FPA2_oriMERGE_base1000_clusterWrapper(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList);
    end

    % FPA for low flux rxns
    changeCobraSolverParams('LP','optTol', 10e-9);
    changeCobraSolverParams('LP','feasTol', 10e-9);
    smallTarget = intersect(lowFluxRxns,targetRxns);
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP_small] = FPA_clusterWrapper(model,smallTarget,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList,...
                                                    {},1,{},0.005);
    elseif strcmp(FPAtype,'new')
        [FP_small] = FPA2_oriMERGE_base1000_clusterWrapper(model,smallTarget,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList,...
                                                    {},1,{},0.005);
    end
      % update small fluxes and wrap up results                             
    for i = 1:size(FP,1)
        if any(strcmp(targetRxns{i},smallTarget))
            ind2 = strcmp(targetRxns{i},smallTarget);
            if FVA_f(strcmp(FVA_label,targetRxns{i})) > 1e-9 && FVA_f(strcmp(FVA_label,targetRxns{i}))<1 
                for j = 1:length(master_expression)
                     FP{i,j}(1) = FP_small{ind2,j}(1);
                end
            end
            if FVA_r(strcmp(FVA_label,targetRxns{i})) < -1e-9 && FVA_r(strcmp(FVA_label,targetRxns{i}))>-1 
                for j = 1:length(master_expression)
                     FP{i,j}(2) = FP_small{ind2,j}(2);
                end
            end
        end
    end
    relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
    relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
    for i = 1:size(FP,1)
        for j = 1:length(master_expression)
            if abs(FP{i,j}(1)) < 1e-9
                FP{i,j}(1) = 0;
            end
            if abs(FP{i,j}(2)) < 1e-9
                FP{i,j}(2) = 0;
            end
            if abs(FP{i,end}(1)) < 1e-9
                FP{i,end}(1) = 0;
            end
            if abs(FP{i,end}(2)) < 1e-9
                FP{i,end}(2) = 0;
            end
            relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
        end
    end
    output_collections.INT = {relFP_f,relFP_r};
	
    % next do the transporters
    FVA_label = model.rxns;
    targetRxns = targetRxns_TSP;
    % the in-organic co-transportable metabilites are:
    helperMet = {'CO2';'AMP';'NADP+';'NADPH';'PPi';'O2';'NADH';'NAD+';
                'Pi';'ADP';'CoA';'ATP';'H2O';'H+';'GTP';'GDP';
                'Electron Transfer Flavoprotein Reduced';'Electron Transfer Flavoprotein Oxidized';
                'L-carnitine';'FAD';'FADH2';'Na+';'hydroxide';'chloride';'iodide'};
            
    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    if strcmp(FPAtype,'original')
        [FP] = MET_FPA_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    elseif strcmp(FPAtype,'new')
        [FP] = MET_FPA2_oriMERGE_base1000_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    end

    % FPA for low flux rxns
    changeCobraSolverParams('LP','optTol', 10e-9);
    changeCobraSolverParams('LP','feasTol', 10e-9);
    smallTarget = intersect(lowFluxRxns,targetRxns);
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP_small] = MET_FPA_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    elseif strcmp(FPAtype,'new')
        [FP_small] = MET_FPA2_oriMERGE_base1000_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    end
    % update small fluxes and wrap up results                             
    for i = 1:size(FP,1)
        if any(strcmp(targetRxns{i},smallTarget))
            ind2 = strcmp(targetRxns{i},smallTarget);
            if FVA_f(strcmp(FVA_label,targetRxns{i})) > 1e-9 && FVA_f(strcmp(FVA_label,targetRxns{i}))<1 
                for j = 1:length(master_expression)
                     FP{i,j}(1) = FP_small{ind2,j}(1);
                end
            end
            if FVA_r(strcmp(FVA_label,targetRxns{i})) < -1e-9 && FVA_r(strcmp(FVA_label,targetRxns{i}))>-1 
                for j = 1:length(master_expression)
                     FP{i,j}(2) = FP_small{ind2,j}(2);
                end
            end
        end
    end
    relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
    relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
    for i = 1:size(FP,1)
        for j = 1:length(master_expression)
            if abs(FP{i,j}(1)) < 1e-9
                FP{i,j}(1) = 0;
            end
            if abs(FP{i,j}(2)) < 1e-9
                FP{i,j}(2) = 0;
            end
            if abs(FP{i,end}(1)) < 1e-9
                FP{i,end}(1) = 0;
            end
            if abs(FP{i,end}(2)) < 1e-9
                FP{i,end}(2) = 0;
            end
            relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
        end
    end
    output_collections.TSP = {relFP_f,relFP_r};
	% finally, do the demands
    load('FVA_cmp.mat');
    load('lowFluxCmp.mat');
    targetRxns = targetRxns_DM;
    FVA_label = targetRxns;
    lowFluxRxns = cellfun(@(x) ['NewMet_',x],lowFluxCmps,'UniformOutput',false);
    % the in-organic co-transportable metabilites are:
    helperMet = {'CO2';'AMP';'NADP+';'NADPH';'PPi';'O2';'NADH';'NAD+';
                'Pi';'ADP';'CoA';'ATP';'H2O';'H+';'GTP';'GDP';
                'Electron Transfer Flavoprotein Reduced';'Electron Transfer Flavoprotein Oxidized';
                'L-carnitine';'FAD';'FADH2';'Na+';'hydroxide';'chloride';'iodide'};
            
    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    if strcmp(FPAtype,'original')
        [FP] = MET_FPA_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    elseif strcmp(FPAtype,'new')
        [FP] = MET_FPA2_oriMERGE_base1000_clusterWrapper(model,targetRxns,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet);
    end

    % FPA for low flux rxns
    changeCobraSolverParams('LP','optTol', 10e-9);
    changeCobraSolverParams('LP','feasTol', 10e-9);
    smallTarget = intersect(lowFluxRxns,targetRxns);
    % Finally, run FPA by simply calling:
    if strcmp(FPAtype,'original')
        [FP_small] = MET_FPA_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    elseif strcmp(FPAtype,'new')
        [FP_small] = MET_FPA2_oriMERGE_base1000_clusterWrapper(model,smallTarget,master_expression,distMat,distMat_raw,labels,n, manualPenalty,manualDist,maxDist,blockList, helperMet,...
                                                    {},1,{},0.001);
    end
    % update small fluxes and wrap up results                             
    for i = 1:size(FP,1)
        if any(strcmp(targetRxns{i},smallTarget))
            ind2 = strcmp(targetRxns{i},smallTarget);
            if FVA_f(strcmp(FVA_label,targetRxns{i})) > 1e-9 && FVA_f(strcmp(FVA_label,targetRxns{i}))<1 
                for j = 1:length(master_expression)
                     FP{i,j}(1) = FP_small{ind2,j}(1);
                end
            end
            if FVA_r(strcmp(FVA_label,targetRxns{i})) < -1e-9 && FVA_r(strcmp(FVA_label,targetRxns{i}))>-1 
                for j = 1:length(master_expression)
                     FP{i,j}(2) = FP_small{ind2,j}(2);
                end
            end
        end
    end
    relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
    relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
    for i = 1:size(FP,1)
        for j = 1:length(master_expression)
            if abs(FP{i,j}(1)) < 1e-9
                FP{i,j}(1) = 0;
            end
            if abs(FP{i,j}(2)) < 1e-9
                FP{i,j}(2) = 0;
            end
            if abs(FP{i,end}(1)) < 1e-9
                FP{i,end}(1) = 0;
            end
            if abs(FP{i,end}(2)) < 1e-9
                FP{i,end}(2) = 0;
            end
            relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
        end
    end
    output_collections.DM = {relFP_f,relFP_r};
end
if ~strcmp(targetType,'HMDB')
    % update small fluxes                              
    for i = 1:size(FP,1)
        if any(strcmp(targetRxns{i},smallTarget))
            ind2 = strcmp(targetRxns{i},smallTarget);
            if FVA_f(strcmp(FVA_label,targetRxns{i})) > 1e-9 && FVA_f(strcmp(FVA_label,targetRxns{i}))<1 
                for j = 1:length(master_expression)
                     FP{i,j}(1) = FP_small{ind2,j}(1);
                end
            end
            if FVA_r(strcmp(FVA_label,targetRxns{i})) < -1e-9 && FVA_r(strcmp(FVA_label,targetRxns{i}))>-1 
                for j = 1:length(master_expression)
                     FP{i,j}(2) = FP_small{ind2,j}(2);
                end
            end
        end
    end
    %% get the final result of reltaive flux potential
    relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
    relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
    for i = 1:size(FP,1)
        for j = 1:length(master_expression)
            if abs(FP{i,j}(1)) < 1e-9
                FP{i,j}(1) = 0;
            end
            if abs(FP{i,j}(2)) < 1e-9
                FP{i,j}(2) = 0;
            end
            if abs(FP{i,end}(1)) < 1e-9
                FP{i,end}(1) = 0;
            end
            if abs(FP{i,end}(2)) < 1e-9
                FP{i,end}(2) = 0;
            end
            relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
        end
    end
    save(['output/FPA_',targetType, '_',expressionDataType,'_',FPAtype,'FPA_',distanceMat,'Dist_order',distanceOrder,'_',networkType,'Network.mat'],'relFP_f','relFP_r');
else
    save(['output/FPA_',targetType, '_',expressionDataType,'_',FPAtype,'FPA_',distanceMat,'Dist_order',distanceOrder,'_',networkType,'Network.mat'],'output_collections');
end
%save('FPA_workspace.mat');

end