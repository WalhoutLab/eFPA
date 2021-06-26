% one significant modification has been made in this version: the low
% category in the intial category is discarded. the a new low category is
% created totally based on DE analysis;
% as a dirty solution for daf2, we use hard foldchange cutoff for DE list,
% so that assume any pairwise DE is very significantly changed. then if
% there is gradiant across four sample, the lowest will go to low, highest
% goes to high, and the two in the middle stay dynamic.
% for a systematic solution of large number integration, a histogram should
% be calculated for each gene, and then the category for this gene in each
% sample shuold be defined by the distribution.
wd = './input/humanModel/GeneCategories/';
%this v2 change the fc calculation to TPM/min(High genes)
%% doing refinement by TS score
TsTbl_RNA = readtable('input/humanModel/RNATSscore.xlsx');
% set NA to 0 (population mean)
tmp = TsTbl_RNA{:,5:end};
tmp(isnan(tmp)) = 0;
TsTbl_RNA{:,5:end} = tmp;
TsTbl_pro = readtable('input/humanModel/proteinTSscore.xlsx');
% set NA to 0 (population mean)
tmp = TsTbl_pro{:,5:end};
tmp(isnan(tmp)) = 0;
TsTbl_pro{:,5:end} = tmp;

TPM = readtable('./input/humanModel/RNATissueMedian_log2TPM.xlsx');
TPM{:,2:end} = 2.^TPM{:,2:end};
load('./input/humanModel/ihuman_COBRA.mat');
% (specific to worm model)
% worm expression data are often labeled with WormBase ID (WBID). So, we
% provided a simple lookup table for iCEL1314.
% lookupTbl = readtable('./input/IDtbl.csv');
 % specifically, we need to convert the gene ID for C. elegans model
% DEtbl = DEtbl(ismember(DEtbl.WBID,lookupTbl.WormBase_Gene_ID),:);%we only keep the genes in the model
% [A,B] = ismember(DEtbl.WBID,lookupTbl.WormBase_Gene_ID);
% DEtbl.GeneID(A) = lookupTbl.ICELgene(B(A));

%% refine the category - for RNA
% when TS score is higher than 2.5 move up; lower than 2.5, move down
names = TPM.Properties.VariableNames(2:end);
TScutoff = 2.5;
for sampleName = names
    sampleName = sampleName{:};
    %% curate the classification
    % make the up/down set
    % the purpose is to fine-tune the dynamic category
    upGenes = intersect(model.genes,TsTbl_RNA.ensembl_id(TsTbl_RNA.(sampleName) >= TScutoff));
    downGenes = intersect(model.genes,TsTbl_RNA.ensembl_id(TsTbl_RNA.(sampleName) <= -TScutoff));
    load([wd,'categ_',sampleName,'.mat']);
    genes = [upGenes;downGenes];
    zero2dyn = 0;
    low2dyn = 0;
    dyn2high = 0;
    dyn2low = 0;
    for i = 1: length(genes)
        mygene  = genes{i};
        if any(strcmp(mygene,upGenes)) 
            if any(strcmp(mygene, ExpCateg.zero)) % zero or low move to dynamic
                fprintf('%s is up but in zero category in %s\n',mygene,sampleName);
                error('check!');
                ExpCateg.zero(strcmp(mygene,ExpCateg.zero)) = [];%detele the old
                ExpCateg.dynamic(end+1,1) = {mygene};
                zero2dyn = zero2dyn +1;
            elseif any(strcmp(mygene, ExpCateg.low)) % zero or low move to dynamic
                fprintf('%s is up but in low category in %s\n',mygene,sampleName);
                ExpCateg.low(strcmp(mygene,ExpCateg.low)) = [];%detele the old
                ExpCateg.dynamic(end+1,1) = {mygene};   
                low2dyn = low2dyn + 1;
            elseif any(strcmp(mygene, ExpCateg.dynamic)) % dynamic to high
                ExpCateg.dynamic(strcmp(mygene,ExpCateg.dynamic)) = [];%detele the old
                ExpCateg.high(end+1) = {mygene};
                dyn2high = dyn2high + 1;
            end   
        else 
            if any(strcmp(mygene, ExpCateg.dynamic)) % dynamic to low
                % this filter could be removed if the network is too
                % similar (or want to be more greedy).
                %if any(strcmp(mygene, ExpCateg.dynamic_low_tail)) % only move genes that expression is not very high
                    ExpCateg.dynamic(strcmp(mygene,ExpCateg.dynamic)) = [];%detele the old
                    ExpCateg.low(end+1,1) = {mygene};   
                    dyn2low = dyn2low + 1;
                %end
            end   
        end
    end
%     ExpCateg = rmfield(ExpCateg,'dynamic_low_tail');
%     ExpCateg = rmfield(ExpCateg,'dynamic_high_tail');
    % save result first
    save([wd,'RNA_refined_',sampleName,'.mat'],'ExpCateg');
    if zero2dyn == 0 &&low2dyn==0&&dyn2high==0&&dyn2low==0
        fprintf('In %s, nothing was changed\n',sampleName);
    else
        if zero2dyn ~= 0
            fprintf('In %s, %d genes were moved from zero to dynamic\n',sampleName, zero2dyn);
        end
        
        if low2dyn ~= 0
        fprintf('In %s, %d genes were moved from low to dynamic\n',sampleName, low2dyn);
        end
        
        if dyn2high~=0
            fprintf('In %s, %d genes were moved from dynamic to high\n',sampleName, dyn2high);
        end
        
        if dyn2low~=0
            fprintf('In %s, %d genes were moved from dynamic to low\n',sampleName, dyn2low);
        end
    end
end
%% check the refinement effectiveness
% check common classification
highComm = model.genes;
dynamicComm = model.genes;
lowComm = model.genes;
zeroComm = model.genes;

highUnion = {};
dynamicUnion = {};
lowUnion = {};
zeroUnion = {};

names =TPM.Properties.VariableNames(2:end);
for i = 1:length(names)
    sampleName = names{i};
    load([wd,'RNA_refined_',sampleName,'.mat'],'ExpCateg'); %change suffix 
    highComm = intersect(ExpCateg.high,highComm);
    dynamicComm = intersect(ExpCateg.dynamic,dynamicComm);
    lowComm = intersect(ExpCateg.low,lowComm);
    zeroComm = intersect(ExpCateg.zero,zeroComm);
    
    highUnion = union(ExpCateg.high,highUnion);
    dynamicUnion = union(ExpCateg.dynamic,dynamicUnion);
    lowUnion = union(ExpCateg.low,lowUnion);
    zeroUnion = union(ExpCateg.zero,zeroUnion);
end
fprintf('common high genes are %d, while union are %d \n',length(highComm),length(highUnion));
fprintf('common dynamics genes are %d, while union are %d \n',length(dynamicComm),length(dynamicUnion));
fprintf('common low genes are %d, while union are %d \n',length(lowComm),length(lowUnion));
fprintf('common zero genes are %d, while union are %d \n',length(zeroComm),length(zeroUnion));

%% refine the category - for protein
exceptionGenes_To_dynamic = {'ENSG00000205186','ENSG00000144035','ENSG00000182591','ENSG00000184210'};% special rules based on manual inspection 

% when TS score is higher than 2.5 move up; lower than 2.5, move down
names = TPM.Properties.VariableNames(2:end);
TScutoff = 2.5;
for sampleName = names
    sampleName = sampleName{:};
    %% curate the classification
    % make the up/down set
    % the purpose is to fine-tune the dynamic category
    upGenes = intersect(model.genes,TsTbl_pro.ensembl_id(TsTbl_pro.(sampleName) >= TScutoff));
    downGenes = intersect(model.genes,TsTbl_pro.ensembl_id(TsTbl_pro.(sampleName) <= -TScutoff));
    load([wd,'categ_',sampleName,'.mat']);
    genes = [upGenes;downGenes];
    zero2dyn = 0;
    low2dyn = 0;
    dyn2high = 0;
    dyn2low = 0;
    for i = 1: length(genes)
        mygene  = genes{i};
        if any(strcmp(mygene,upGenes)) 
            if any(strcmp(mygene, ExpCateg.zero)) % zero or low move to dynamic
                if any(strcmp(exceptionGenes_To_dynamic,mygene))
                    ExpCateg.zero(strcmp(mygene,ExpCateg.zero)) = [];%detele the old
                    ExpCateg.dynamic(end+1,1) = {mygene};
                    zero2dyn = zero2dyn +1;
                else
                    fprintf('%s is up but in zero category in %s\n',mygene,sampleName);
                    error('check!');
                end
            elseif any(strcmp(mygene, ExpCateg.low)) % zero or low move to dynamic
                fprintf('%s is up but in low category in %s\n',mygene,sampleName);
                ExpCateg.low(strcmp(mygene,ExpCateg.low)) = [];%detele the old
                ExpCateg.dynamic(end+1,1) = {mygene};   
                low2dyn = low2dyn + 1;
            elseif any(strcmp(mygene, ExpCateg.dynamic)) % dynamic to high
                ExpCateg.dynamic(strcmp(mygene,ExpCateg.dynamic)) = [];%detele the old
                ExpCateg.high(end+1) = {mygene};
                dyn2high = dyn2high + 1;
            end   
        else 
            if any(strcmp(mygene, ExpCateg.dynamic)) % dynamic to low
                % this filter could be removed if the network is too
                % similar (or want to be more greedy).
                %if any(strcmp(mygene, ExpCateg.dynamic_low_tail)) % only move genes that expression is not very high
                    ExpCateg.dynamic(strcmp(mygene,ExpCateg.dynamic)) = [];%detele the old
                    ExpCateg.low(end+1,1) = {mygene};   
                    dyn2low = dyn2low + 1;
                %end
            end   
        end
    end
%     ExpCateg = rmfield(ExpCateg,'dynamic_low_tail');
%     ExpCateg = rmfield(ExpCateg,'dynamic_high_tail');
    % save result first
    save([wd,'protein_refined_',sampleName,'.mat'],'ExpCateg');
    if zero2dyn == 0 &&low2dyn==0&&dyn2high==0&&dyn2low==0
        fprintf('In %s, nothing was changed\n',sampleName);
    else
        if zero2dyn ~= 0
            fprintf('In %s, %d genes were moved from zero to dynamic\n',sampleName, zero2dyn);
        end
        
        if low2dyn ~= 0
        fprintf('In %s, %d genes were moved from low to dynamic\n',sampleName, low2dyn);
        end
        
        if dyn2high~=0
            fprintf('In %s, %d genes were moved from dynamic to high\n',sampleName, dyn2high);
        end
        
        if dyn2low~=0
            fprintf('In %s, %d genes were moved from dynamic to low\n',sampleName, dyn2low);
        end
    end
end
%% check the refinement effectiveness
% check common classification
highComm = model.genes;
dynamicComm = model.genes;
lowComm = model.genes;
zeroComm = model.genes;

highUnion = {};
dynamicUnion = {};
lowUnion = {};
zeroUnion = {};

names =TPM.Properties.VariableNames(2:end);
for i = 1:length(names)
    sampleName = names{i};
    load([wd,'protein_refined_',sampleName,'.mat'],'ExpCateg'); %change suffix 
    highComm = intersect(ExpCateg.high,highComm);
    dynamicComm = intersect(ExpCateg.dynamic,dynamicComm);
    lowComm = intersect(ExpCateg.low,lowComm);
    zeroComm = intersect(ExpCateg.zero,zeroComm);
    
    highUnion = union(ExpCateg.high,highUnion);
    dynamicUnion = union(ExpCateg.dynamic,dynamicUnion);
    lowUnion = union(ExpCateg.low,lowUnion);
    zeroUnion = union(ExpCateg.zero,zeroUnion);
end
fprintf('common high genes are %d, while union are %d \n',length(highComm),length(highUnion));
fprintf('common dynamics genes are %d, while union are %d \n',length(dynamicComm),length(dynamicUnion));
fprintf('common low genes are %d, while union are %d \n',length(lowComm),length(lowUnion));
fprintf('common zero genes are %d, while union are %d \n',length(zeroComm),length(zeroUnion));

