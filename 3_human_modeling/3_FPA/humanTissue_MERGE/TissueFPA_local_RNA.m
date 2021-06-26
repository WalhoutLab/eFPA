
distanceMat = 'weighted';
FPAtype = 'new';
expressionDataType = 'RNA_TS_common';
networkType = 'tissue';
targetType = 'rxn';
n = 6;
targetRxns = {'HMR_2182'}; % GABA

%% This is to reproduce the C. elegans tissue FPA results in supp table S7 and S8
%% set up the env variables
addpath ~/cobratoolbox/
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
initCobraToolbox(false)
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
% Users can uncomment the following command to load from the distance 
% calculator output. Here we load directly from saved matlab variable,
% because of the file size restriction of GitHub
% distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix_recon2_2.txt','FileType','text','ReadRowNames',true); % we load from the output of the distance calculator. For usage of distance calculator, please refer to the MetabolicDistance folder
if strcmp(distanceMat,'original')
    distance_raw = readtable('./input/distanceMatrix.txt','FileType','text','ReadRowNames',true); % we load from the output of the distance calculator. For usage of distance calculator, please refer to the MetabolicDistance folder
elseif strcmp(distanceMat,'weighted')
    distance_raw = readtable('./input/distanceMatrix_weighted.txt','FileType','text','ReadRowNames',true); % we load from the output of the distance calculator. For usage of distance calculator, please refer to the MetabolicDistance folder
    % there are some NaN. seems related to transport rxns; unknown why. assume
    % inf for now! 11252020
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

% ```set the special penalties (if desired)```
extRxns = model.rxns(findExcRxns_XL(model));
% sideUptake = model.rxns(model.S(strcmp(model.mets,'sideNutr'),:)~=0);
% majorUptake = model.rxns(model.S(strcmp(model.mets,'majorNutr'),:)~=0);

manualPenalty = [extRxns];%;sideUptake;majorUptake];
manualPenalty(:,2) = [mat2cell(zeros(length(extRxns),1),ones(length(extRxns),1))];%;...
                      %mat2cell(10.*ones(length(sideUptake),1),ones(length(sideUptake),1));...
                      %mat2cell(zeros(length(majorUptake),1),ones(length(majorUptake),1))];
                  
% special distance for side - to penalize distal side usage
manualDist = {};%sideUptake;
%manualDist(:,2) = mat2cell(10.*ones(length(sideUptake),1),ones(length(sideUptake),1));
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
parpool(4,'SpmdEnabled',false); % we run this code on a 40-core lab server
%% PART I: FPA for normal reactions (reaction-centered FPA)
penalty = calculatePenalty(model,master_expression,manualPenalty);
%% PART I: FPA for normal reactions (reaction-centered FPA)
if strcmp(targetType,'rxn')
    load('lowFluxRxns.mat');
    changeCobraSolverParams('LP','optTol', 10e-7);
    changeCobraSolverParams('LP','feasTol', 10e-7);
    % FPA with normal alpha
    % Finally, run FPA by simply calling:
    if ~strcmp(targetRxns,lowFluxRxns)
        if strcmp(FPAtype,'original')
            [FP,FP_solutions] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList, {},0,penalty);
        elseif strcmp(FPAtype,'new')
            [FP,FP_solutions] = FPA2_oriMERGE_base2(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,manualDist,maxDist,blockList, {},0,penalty);
        end
    else
        error('low flux rxn');
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

%% check rFP
figure(1)
c = categorical(conditions);
bar(c, relFP_f(1,:))
ylabel('raw')
ylim([0,1])
figure(2)
c = categorical(conditions);
bar(c, normalize(relFP_f(1,:),'center','median'))
ylabel('raw')
ylim([-0.5,1])
printGPRForRxns(model,targetRxns);
%%
figure(1)
c = categorical(conditions);
bar(c, relFP_r(1,:))
ylabel('raw')
ylim([0,1])
figure(2)
c = categorical(conditions);
bar(c, normalize(relFP_r(1,:),'center','median'))
ylabel('raw')
ylim([-0.5,1])
printGPRForRxns(model,targetRxns);
%% check expression 
name2 = 'ENSG00000151093';
%expTbl(strcmp(expTbl.ensembl_id,name2),:)
figure(3)
c = categorical(conditions);
bar(c, expTbl{strcmp(expTbl.ensembl_id,name2),5:end})
ylabel('expression level')
%% check a genes' impact
model.rxns(logical(model.rxnGeneMat(:,strcmp(model.genes,name2))))
%% check distribution
model_irr = convertToIrreversible(model);
%%
tissue = 'Liver';
FP_solutions{1,strcmp(conditions,tissue)}{1}.rxns = model_irr.rxns;
mytbl = listFPtbl(model_irr,FP_solutions{1,strcmp(conditions,tissue)}{1});
%%
tissue = 'MuscleSkeletal';
FP_solutions{1,strcmp(conditions,tissue)}{2}.rxns = model_irr.rxns;
mytbl = listFPtbl(model_irr,FP_solutions{1,strcmp(conditions,tissue)}{2});
%% track flux distirbution 
mytbl2 = listRxn(model_irr,FP_solutions{1,strcmp(conditions,tissue)}{1}.x,'pyruvate [Extracellular]');

%%
myrxn = {'HMR_2181'};
figure(7)
c = categorical([conditions,{'super'}]);
bar(c, penalty(strcmp(model.rxns,myrxn{:}),1:end))
ylabel('rxn-level penalty')
%% some notes
% tyr -> dopa -> dopamine is captured 
% tyr -> dopa --> dopamin --> noradrenaline is NOT captured. dopamine is
% taken up directly 
% --> 5hp --> serotonin is captured but the expression of the first step is mostly not
% measured 

% 'HMR_2259' to 'HMR_2265' 'Fatty acid elongation (odd-chain)' is a
% pathway-level prediction 