%% This walkthrough guidance will take user through the application of FPA on a generic model of C. elegans, or other metabolic models (we used human model recon2.2 as an example).
%
% Please refer to: 
% Yilmaz, Li et al., 2020, in review.

%% CAUTION ABOUT MATLAB VERSION
% This demo is developed and tested in MATLAB R2019a version. An earlier
% version may encounter errors with "readtable" and "writematrix"
% functions.

%% Add required path 
% **** make sure you do this every time you run the demo!****
% add path for required functions/inputs
addpath ./input/
addpath ./scripts/
addpath ./../input/
addpath ./../bins/
initCobraToolbox(false);


%% Part I: APPLICATION TO THE GENERIC C. ELEGANS MODEL
%% 1. Prepare the model
load('input/ihuman_COBRA.mat');
mkdir('distance_inputs');
%% 2. (OPTIONAL) generate the distance matrix (skip this if you already have the matrix)
% This section will guide users to generate the input for metabolic
% distance calculator from a COBRA model.
writematrix(model.S,'distance_inputs/Smatrix_regular.txt');
writecell(model.rxns,'distance_inputs/reactions_regular.txt');
writecell(model.mets,'distance_inputs/metabolites_regular.txt');
writematrix(model.lb,'distance_inputs/LB_regular.txt');
writematrix(model.ub,'distance_inputs/UB_regular.txt');
% Users need to determine their own byproduct metabolites (for details, see Structured Methods in our paper) 
% byProducts = {'co2';'amp';'nadp';'nadph';'ppi';'o2';'nadh';'nad';'pi';'adp';'coa';'atp';'h2o';'h';'gtp';'gdp';'etfrd';'etfox';'crn';'fad';'fadh2'};% the typical set of byproducts
byProducts = {'CO2';'AMP';'NADP+';'NADPH';'PPi';'O2';'NADH';'NAD+';
            'Pi';'ADP';'CoA';'ATP';'H2O';'H+';'GTP';'GDP';
            'Electron Transfer Flavoprotein Reduced';'Electron Transfer Flavoprotein Oxidized';
            'L-carnitine';'FAD';'FADH2';'Na+'};% the typical set of byproducts
% Add compartment label to byproducts
byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,' \[.+\]$',''),model.metNames, 'UniformOutput',false),byProducts));
writecell(byProducts,'distance_inputs/byproducts_regular.txt');
% Then you can use the output files in folder the "distance_inputs" folder 
% for calculating distances. Please follow the Distance calculator section 
% in Github
%% For an earlier MATLAB version, user may consider the following codes:
% dlmwrite('distance_inputs/Smatrix_regular.txt',full(model.S),',');
% fid=fopen('distance_inputs/reactions_regular.txt','w');
% fprintf(fid,'%s\n',model.rxns{:});
% fclose(fid);
% fid=fopen('distance_inputs/metabolites_regular.txt','w');
% fprintf(fid,'%s\n',model.mets{:});
% fclose(fid);
% dlmwrite('distance_inputs/LB_regular.txt',model.lb,',');
% dlmwrite('distance_inputs/UB_regular.txt',model.ub,',');
% fid=fopen('distance_inputs/byproducts_regular.txt','w');
% fprintf(fid,'%s\n',byProducts{:});
% fclose(fid);

%% 3. Load the expression files, distance matrix, and other optional inputs

% ```load expression files```
% Expression matrix can be in plain text and in any normalized 
% quantification metric like TPM or FPKM.
expTbl = readtable('./../1_IMAT++/input/humanModel/NX/rna_tissue_consensus_metabolic.csv');
% This table provides normalized expression (NX) values for metabolic genes 
% across human tissues. It is extracted from "RNA consensus tissue gene 
% data" downloaded from www.proteinatlas.org as of Jun-1-2020

% ```preprocess the expression table```
% Since expression tables from different data source may have different
% format, we require user to re-organize your table according to the
% following procedure.
%
% In FPA, we convert the expression table into a cell array of 
% structure variables called "master_expression".
% REQUIREMENTS OF "master_expression"
% (1) it must be a cell arrary of structure variables;
% (2) each structure variable contains the expression information of one
%     condition; and the order of the structure variables in the cell array 
%     determines the order of Flux Potential (FP) values in the output;
% (3) each structure variable must contain two fields: "genes" amd "value".
%     "genes" and "value" should be of equal length.
% (4) each of the input "genes" should be measured in ALL conditions. In 
%     other words, it is not recommended to have different "genes" list in
%     different conditions. 
%
% Make the master_expression for human tissues
master_expression = {};
conditions = unique(expTbl.Tissue); % conditions to calculate FPA on
% Determine the genes that are measured in all conditions (requirement #4)
targetGenes = intersect(expTbl.HGNC(strcmp(expTbl.Tissue,conditions{1})), model.genes); % genes both detected and in the model
for i = 2: length(conditions)
    targetGenes = intersect(expTbl.HGNC(strcmp(expTbl.Tissue,conditions{i})),targetGenes);
end
% Fill in the "master_expression" cell array
for i = 1:length(conditions) 
    expression = struct();
    expression.genes = targetGenes; 
    tissueInd = strcmp(expTbl.Tissue,conditions{i});
    geneLabel = expTbl.HGNC(tissueInd);
    NXval = expTbl.NX(tissueInd);
    % reorder
    [A B] = ismember(targetGenes,geneLabel);
    expression.value = NXval(B(A));
    master_expression{i} = expression;
end


% ```load the distance matrix```
% Users can uncomment the following command to load from the distance 
% calculator output. Here we load directly from saved matlab variable,
% because of the file size restriction of GitHub
% distance_raw = readtable('./../MetabolicDistance/Output/distanceMatrix_recon2_2.txt','FileType','text','ReadRowNames',true); % we load from the output of the distance calculator. For usage of distance calculator, please refer to the MetabolicDistance folder
load('./input/distance_raw_recon2_2.mat');
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
% *** this could be skipped for quick and simple test of FPA ***
% FPA allows user to avoid the penalty for desired reactions. This could
% be done by setting special penalties (`manualPenalty` parameter). 
%
% In general, we recommend to avoid the penalty of external reactions in FPA.
% Therefore, we set penalty for all `Exchange`, `Demand`, and `Sink` reactions to 0. 
% Users may need to interactively tune their special penalties for the best flux
% distribution in FPA calculation
extRxns = model.rxns(cellfun(@(x) ~isempty(regexp(x,'^(EX_|sink_|DM_)', 'once')),model.rxns));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

%% 4. run regular FPA analysis
% As part of the MERGE package, we recommend user to integrate the result
% of iMAT++ to the FPA analysis. That's saying, to block all reactions that
% don't carry flux in the feasible solution space. These reactions are
% identified by FVA analysis conjoined with IMAT++. Please refer to 
% walkthrough tutorial of IMAT++ and "1_IMAT++/walkthrough_large_scale_FVA.m" 
% for getting the FVA result.
%
% Assuming the FVA is done, we should have the level table ready for each 
% condition (we generated the level tables in the FVA walkthrough). For
% demo purpose, we only performed FVA for 17 selected tissues. Therefore,
% we will only perform FPA on these 17 tissues.

% ```setup `master_expression` for the 17 selected tissues```
ExampleTissues = {'cerebral cortex','spinal cord','midbrain','pons and medulla',...%neuronal tissues (<--> worm neuron/glia)
                'duodenum','colon','small intestine',... % digestive tissues (<--> worm intestine)
                'ovary','seminal vesicle','testis',... % reproductive tissues (<--> worm gonad)
                'heart muscle','smooth muscle','tongue','skeletal muscle',... (<--> worm muscle, pharnyx)
                'skin','kidney','liver'};  %(skin and liver <--> worm hypodermis)
[A B] = ismember(ExampleTissues,conditions);
master_expression = master_expression(B(A));
% FPA will be reported in the same order as that of the tissues in 
% the `ExampleTissues`

% ```setup the block list (context-specific networks)```
% First, let's load and merge the level tables for each tissue
ExampleTissues_name2 = {'cerebralCortex','spinalCord','midbrain','ponsAndMedulla',...%neuronal tissues (<--> worm neuron/glia)
                'duodenum','colon','smallIntestine',... % digestive tissues (<--> worm intestine)
                'ovary','seminalVesicle','testis',... % reproductive tissues (<--> worm gonad)
                'heartMuscle','smoothMuscle','tongue','skeletalMuscle',... (<--> worm muscle, pharnyx)
                'skin','kidney','liver'};  %(skin and liver <--> worm hypodermis)
for i = 1:length(ExampleTissues_name2)
    load(['./../1_iMAT++/output/humanTissue/FVA_levels/',ExampleTissues_name2{i},'_levels.mat']);
    levelTbl_f(:,i) = levels_f;
    levelTbl_r(:,i) = levels_r;
end
% Then, get the list of to-block reactions
% because the FPA is done using the irreversible model, and therefore, 
% the rxnID is different from the original. We provided a function to get 
% the new IDs.
blockList = getBlockList(model,levelTbl_f,levelTbl_r);
%
% IMPORTANT: the order of tissues (or conditions) in `master_expression` 
% must be the SAME as in `blockList`. Please notice that we aligned the 
% order when we load the FVA outputs above. 


% ```setup some basic parameters for FPA```
n = 1.5; % distance order
changeCobraSolverParams('LP','optTol', 10e-9);
changeCobraSolverParams('LP','feasTol', 10e-9);

%% we perform FPA analysis for four targets as a demo
% (1) isocitrate dehydrogenase reaction (reaction-centric)
% (2) alpha-ketoglutarate dehydrogenase reaction (reaction-centric)
% (3) the cellular demand of melanin (metabolite-centric)
% (4) GABA transporter (metabolite-centric)

%% REACTION-LEVEL ANALYSIS
%% (1) and (2) 
% Reaction-level analysis is simple. 
% Set the target reactions for FPA:
targetRxns = {'ICDHyrm','AKGDm'};
% Finally, run FPA by simply calling:
[FP,FP_solutions] = FPA(model,targetRxns,master_expression,distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);

%% METABOLITE-LEVEL ANALYSIS
% Metabolite-level analysis can be performed in two ways: 
% (i) on a metabolite sink/demand reaction or 
% (ii) on a metabolite transporter reaction.
% We illustrate these individually.

%% (3) the cellular demand of melanin (on a metabolite demand reaction)
% First, let's create the new demand reaction for metabolite-level analysis of melanin
% (You don't need to add the new reaction if the target demand/sink already exists in the model)
myMet = 'melanin[c]';
myRxnName = ['DMN',myMet];
model_new = addReaction(model,myRxnName,'reactionFormula',[myMet,' -->'],'geneRule', 'NA','printLevel',0);
%
% Then, we need to update the distance matrix because of adding the new demand reaction
[distMat_raw_new,labels_new] = updateDistance(myRxnName,model_new,distMat_raw,labels,max(max(distMat_raw)));   
distMat_new = distMat_raw_new;
for ii = 1:size(distMat_new,1)
    for jj = 1:size(distMat_new,2)
        % take the minimal of the raw distance (see methods for details)
        distMat_new(ii,jj) = min([distMat_raw_new(ii,jj),distMat_raw_new(jj,ii)]);
    end
end
%
% Finally, call FPA function again
targetRxns = {myRxnName};
[FP(3,:),FP_solutions(3,:)] = FPA(model_new,targetRxns,master_expression,distMat_new,labels_new,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);
%
% Please note, you may need to block the reactions that could uptake your
% target metabolite (i.e, some exchange reaction). These reactions cause 
% shortcut flux in FPA (see methods for explanation). In this melanin 
% example, we don't need to block since there is no such reaction.
% We will illustrate blocking shortcuts in the next example.
 
%% (4) GABA transporter (on a metabolite transporter reaction)
% We use the canonical sodium-dependent GABA transporter as the target transporter
targetRxns = {'ABUTt4_2_r'}; 
% block shortcuts
shortcutRxns = {'ABUTt2r','GABAVESSEC','GABABGTtc','DM_4abut_n_'}; % other GABA transporters and demand reactions that could interfere with the analysis of the target transporter
model_block = changeRxnBounds(model,shortcutRxns,0,'b');
% To do the analysis of GABA degradation potential, we allow the GABA exchange
model_block = changeRxnBounds(model_block,'EX_4abut_e_',-1000,'l');
%
% Finally, call FPA
[FP(4,:),FP_solutions(4,:)] = FPA(model_block,targetRxns,master_expression,distMat,labels,n, manualPenalty,{},max(distMat(~isinf(distMat))),blockList);
%
% Please be advised that you can also analyze the import/export potential
% of other GABA transporters (which we blocked in this case). We only show
% the analysis of one GABA transporter (ABUTt4_2_r) as an example.

%% 5. Calculate relative flux potential (rFP)
relFP_f = nan(size(FP,1),length(master_expression));% flux potential for forward rxns
relFP_r = nan(size(FP,1),length(master_expression));% flux potential for reverse rxns
for i = 1:size(FP,1)
    for j = 1:length(master_expression)
        relFP_f(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
        relFP_r(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
    end
end
%
% The relative flux potentials (rFP) are stored in two variables, relFP_f 
% and relFP_r. '_f' and '_r' represent forward direction and reverse 
% direction, respectively. Rows of the rFP variables are the four queried 
% reactions; columns are the 17 queried tissues. For example, relFP_f(1,1) 
% is the rFP of the forward direction of ICDHym, for the first tissue, 
% cerebral cortex.

%% IMPORTANT NOTICE ON METABOLITE-LEVEL FPA
% If you are doing large-scale metabolite-level FPA analysis,
% please refer to the corresponding section in "TissueFPA.m" on how to
% programmatically block shortcuts and analyze transporter flux potentials.
% (also on how to avoid mapping expression to penalties repeatedly)
% This demo runs metabolite-level FPA in a simplified way, so that the shortcuts 
% are blocked manually. 
% Blocking shortcuts programmatically is a model-specific task, as
% different models have different nomenclature system and symbols. 
% Therefore, we can't provide a universal function to do the job. Users may 
% want to engineer our code for worm tissue FPA to do large-scale
% metabolite-level FPA on their own model.

%% (OPTIONAL) run FPA on the full network
% Building context-sepcific network is costly and sometimes the user may
% want to give a quick trial of FPA. Therefore, it is also an option to run
% FPA as a standalone algorithm without the "blocklist" input. Please refer
% to the Appendix Supplementary Methods section of our paper for pros and 
% cons of FPA on the full network.
%
% To run FPA on the full network, simply write:
% example (1) and (2)
[FP_naive,FP_solutions_naive] = FPA(model,{'ICDHyrm','AKGDm'},master_expression,distMat,labels,n, manualPenalty);
% example (3)
[FP_naive(3,:),FP_solutions_naive(3,:)] = FPA(model_new,{'DMNmelanin[c]'},master_expression,distMat_new,labels_new,n, manualPenalty);
% example (4)
[FP_naive(4,:),FP_solutions_naive(4,:)] = FPA(model_block,{'ABUTt4_2_r'},master_expression,distMat,labels,n, manualPenalty);
% Calculate relative potential
relFP_f_naive = nan(size(FP_naive,1),length(master_expression));% flux potential for forward rxns
relFP_r_naive = nan(size(FP_naive,1),length(master_expression));% flux potential for reverse rxns
for i = 1:size(FP_naive,1)
    for j = 1:length(master_expression)
        relFP_f_naive(i,j) = FP_naive{i,j}(1) ./ FP_naive{i,end}(1);
        relFP_r_naive(i,j) = FP_naive{i,j}(2) ./ FP_naive{i,end}(2);
    end
end

%% 6. Interpretation of FPA results
figure(1)
c = categorical(ExampleTissues);
bar(c,relFP_f(1,:))
title('Isocitrate dehydrogenase (NADP+)')

figure(2)
c = categorical(ExampleTissues);
bar(c,relFP_f(2,:))
title('alpha-ketoglutarate dehydrogenase (NADP+)')

figure(3)
c = categorical(ExampleTissues);
bar(c,relFP_f(3,:))
title('Melanin production')

figure(4)
c = categorical(ExampleTissues);
bar(c,relFP_r(4,:))
title('GABA export')

figure(5)
c = categorical(ExampleTissues);
bar(c,relFP_f(4,:))
title('GABA import')

% INTERPRETATION:

% (1) rFP predictions
% From the rFP value, we see the superiority of skeletal muscle (and other
% muscles) in TCA cycle flux potential (via ICDHym and AKGDm), and skin in
% producing melanin. 
% Similarly, the import and export potential of GABA are also in line with 
% our understanding of human physiology. Importantly, FPA could distinguish 
% the export and import potential of GABA in different tissues. For example,
% neuronal tissues are better at exporting (producing) GABA than liver
% and kidney. However, their importing (degrading) potentials are similar
% to that of liver and kidney. This is in line with manual inspection of
% gene expression data. The GABA synthesis gene, GAD1, is enriched
% specifically in brain; However, the ABAT gene, which can both synthesize
% and degrade GABA, is enriched in both liver, kidney and brain. 

% Finally, by inspecting the flux distribution of ICDHyrm FPA (see NOTE#1), 
% we find that FPA successfully integrates the local expression information 
% of TCA cycle (with a shunt) in analyzing the flux potential of isocitrate 
% dehydrogenase. The flux route is pyr --> accoa --> cit --> icit --> akg 
% --> mal --> oaa --> cit. Thus, FPA successfully made a pathway-level 
% prediction.


% (2) The effect of using block list (context-specific networks)
% We can find that the tissue-specific networks blocks the flux potential 
% of melanin production in all tissues except for skin. This provides 
% clearer prediction that melanin is only produced in skin. More 
% importantly, it increases the resolution of GABA import/export potential
% in different tissues. By kicking out tissues whose network doesn't
% support the production or degredation of GABA, we can distinguish low
% production/consumption from no production/consumption. Besides, we can 
% exlude tissues who express GABA transporter but cannot produce or consume
% GABA metabolically.

%% TECHNICAL NOTES FOR USING FPA

%% 1. To inspect the flux distribution for reported FP values
% the flux distribution is reported for the irreversible model, and is 
% stored in the "FP_solutions". 
% To inspect the flux distribution, first get irreversible model
model_irrev = convertToIrreversible(model);
% Then, we provide a flux tracker for easy-inspection.
mytbl = listRxn(model_irrev,FP_solutions{1,11}{1}.full,'icit[m]');
% mytbl: column#1: rxnID, col#2: rxn flux, col#3: rxn formula, col#4: flux
% contribution to the queried metabolite.
% You can view the "mytbl" variable for flux going in and out of isocitrate 
% in the FPA calculation of heart muscle for ICDHym. 

%% 2. notice on gene names and expression data
% (1) Gene Name Rules:
%   (a) we allow letters, numbers, dot, dash and colon in gene names. Any  
%       other symbol needs to be added in the regexp function of line 43 in 
%       eval_gpr.m.
%   (b) if a gene name appears multiple times in the expression table, we 
%       use the sumation of all expression levels in the GPR parsing step. 

% (2) Expression Data Type:
% By default, we assume that the RNA-seq profile is used as expression
% input. Therefore, we, by default, add a pseudocount of 1 to all the genes,  
% to offset the noise when a gene is lowly expressed. If one would like to 
% use other pseudocount value, or not add pseudocount, please modify
% ./scripts/calculatePenalty.m line 48. Please be advised that we DO NOT
% recommend adding pseudocount for microarray dataset!

% (3) If a gene is undetected, this gene is ignored in the GPR parsing
% step, but other detected genes associated with the same reaction will  
% be still used to calculate the reaction-level penalty.
