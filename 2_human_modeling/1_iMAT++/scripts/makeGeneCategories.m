%% Overview
% This is an interactive gene category builder that helps user to fit and
% build the categories for iMAT++. The expression quantifications (TPM
% or FPKM) could be provided as a csv table.
% it is a equivalent of CatExp algorithm in Yilmaz at al 2020. We made
% minor modifications here. 
% ..Author: Xuhang Li, Mar 2020
mkdir(['input/humanModel/GeneCategories']);
%% part I: load the expression data and adjust the format
TPM = readtable('./input/humanModel/RNATissueMedian_log2TPM.xlsx');
TPM{:,2:end} = 2.^TPM{:,2:end};
load('./input/humanModel/ihuman_COBRA.mat');

metGenesInd = ismember(TPM.gene_id,model.genes);
geneInfo = readtable('./input/humanModel/humanGeneInfo.xlsx');
proteinCodingGenes = geneInfo.EnsemblGeneID(strcmp(geneInfo.LocusType,'gene with protein product'));
proteinCodingInd = ismember(TPM.gene_id,proteinCodingGenes);
%% part II: the classification based on absolute value
fitData = TPM{proteinCodingInd,2:end};
fitData = fitData(:);
%fitData = fitData(:);
histogram(log2(fitData))
% fit bimodel guassian
rng(1126)
x = log2(fitData);% this automatically ignored all the zeros
x(isinf(x)) = [];
fit = fitgmdist(x,3,'Options',statset('Display','final','MaxIter',3000,'TolFun',1e-9),'Replicates',10);
% note: the sigma in output is sigma^2
%% visualization of the guassian distribution fitted
% user can inspect how well the fitting is and whether it seperates into
% two subpopulation as expected.
bins = -15:.5:15;
h = bar(bins,histc(x,bins)/(length(x)*.5),'histc');
h.FaceColor = [.9 .9 .9];
xgrid = linspace(1.1*min(x),1.1*max(x),200);
pdfgrid = pdf(fit,xgrid');
hold on
plot(xgrid,pdfgrid,'-')
hold off
xlabel('x')
ylabel('Probability Density')
%% label the desired cutoff
xline(fit.mu(3),'--r');
fit.mu(3);
xline(fit.mu(2),'--r');
fit.mu(2);
xline(fit.mu(1),'--r');
fit.mu(1);

xline(fit.mu(3) + 2*sqrt(fit.Sigma(3)),'--k');

% though we find that the mean of the two subpopulation gives best
% categories in the datasets we tested, we recommend users to interactively
% evaluate their dataset and find the best thresholds
%% build the gene catagories
zero2low = min(fit.mu) %set thresholds
low2dynamic = median(fit.mu) %set thresholds
dynamic2high = max(fit.mu) %set thresholds
names = TPM.Properties.VariableNames(2:end);%choose all the sample names
metgenes = model.genes;
TPM = TPM(ismember(TPM.gene_id,model.genes),:);%we only keep the genes in the model
for myName = names
    myTPM = log2(TPM.(myName{:}));
    GeneID = TPM.gene_id;
    ExpCateg.zero = GeneID(myTPM < zero2low);
    ExpCateg.low = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
    ExpCateg.dynamic = GeneID(myTPM >= low2dynamic & myTPM < dynamic2high);
    ExpCateg.high = GeneID(myTPM >= dynamic2high);
    % the uncalled genes (i.e., NA and ND)are in dynamic (moderately expressed)
    ExpCateg.dynamic = [ExpCateg.dynamic; metgenes(~ismember(metgenes,GeneID))];
    save(['input/humanModel/GeneCategories/categ_',myName{:},'.mat'],'ExpCateg');
end
%% analyze if the thereshold makes sense
% A critical signature of good gene categories is that most metabolic genes are
% called as high when merging the high categories for all conditions
% (assuming the metabolism across conditions varies), while only a core set
% of genes are high in all of the conditions. This indicates that the
% categorization captures the expression regulation of metabolic genes
% across conditions. Similarly, we expect the zero category to have the
% similar feature. 
myTPM = log2(TPM.(names{1}));
GeneID = TPM.gene_id;
ZeroInAll = GeneID(myTPM < zero2low);
ZeroMerge = GeneID(myTPM < zero2low);
LowInAll = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
LowMerge = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
HighInAll = GeneID(myTPM >= dynamic2high);
HighMerge = GeneID(myTPM >= dynamic2high);
N_zero = [];
N_low = [];
N_dynamic = [];
N_high = [];
for myName = names
    myTPM = log2(TPM.(myName{:}));
    GeneID = TPM.gene_id;
    ExpCateg.zero = GeneID(myTPM < zero2low);
    ExpCateg.low = GeneID(myTPM >= zero2low & myTPM < low2dynamic);
    ExpCateg.dynamic = GeneID(myTPM >= low2dynamic & myTPM < dynamic2high);
    ExpCateg.high = GeneID(myTPM >= dynamic2high);
    % the uncalled genes (i.e., NA and ND)are in dynamic (moderately expressed)
    ExpCateg.dynamic = [ExpCateg.dynamic; metgenes(~ismember(metgenes,GeneID))];
    ZeroInAll = intersect(ZeroInAll,ExpCateg.zero);
    ZeroMerge = union(ZeroMerge,ExpCateg.zero);
    LowInAll = intersect(LowInAll,ExpCateg.low);
    LowMerge = union(LowMerge,ExpCateg.low);
    HighInAll = intersect(HighInAll,ExpCateg.high);
    HighMerge = union(HighMerge,ExpCateg.high);
    N_zero = [N_zero;length(ExpCateg.zero)];
    N_low = [N_low;length(ExpCateg.low)];
    N_dynamic = [N_dynamic;length(ExpCateg.dynamic)];
    N_high = [N_high;length(ExpCateg.high)];
end
fprintf('%d/%d are highly expressed genes in all conditions\n',length(HighInAll),length(model.genes));
fprintf('%d/%d are highly expressed genes in at least one condition\n',length(HighMerge),length(model.genes));
fprintf('%d/%d are lowly expressed genes in all conditions\n',length(LowInAll),length(model.genes));
fprintf('%d/%d are lowly expressed genes in at least one conditions\n',length(LowMerge),length(model.genes));
fprintf('%d/%d are rarely expressed genes in all conditions\n',length(ZeroInAll),length(model.genes));
fprintf('%d/%d are rarely expressed genes in at least one conditions\n',length(ZeroMerge),length(model.genes));
% we offer an additional QC figure for category making (similar to fig. 1E
figure(3)
stackN = [N_zero,N_low,N_dynamic,N_high];
bar(1:size(stackN,1),stackN,'stacked')
xlabel('condition No.')
legend({'zero','low','dynamic','high'});