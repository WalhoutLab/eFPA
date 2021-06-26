%% notes
% flux is in mmole/hr/gDW by chemostat info: gDCW/ml. 
% relative flux is by normalize per cell flux by the dilution rate 
%% load the model
addpath('./../scripts/')
model = loadYeatModel();
%% load the expression files, distance matrix, and other optional inputs
% load proteomics
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
% we use the non-log level
proTbl{:,2:end} = 2.^proTbl{:,2:end};
% normalizaed
coeffTbl = readtable('./../input/YeastJoshua/originalDataTbl/coefficients_gram_per_gDW.xlsx');
% preprocess the expression table
% to facilate the future use of the expression of many samples, we
% re-organize it into a structure variable.
% the FPA matrix will be in the same order as the master_expression
conditions = proTbl.Properties.VariableNames(2:end);
% make a new master_expression for these four conditions.
master_expression_perPro = {};% we call this variable "master_expression"
master_expression_perDW= {};% we call this variable "master_expression"
geneInd = ismember(proTbl.Gene, model.genes); % get the index of genes in the model
for i = 1:length(conditions)
    expression = struct();
    expression.genes = proTbl.Gene(geneInd);
    expression.value = proTbl.(conditions{i})(geneInd);
    master_expression_perPro{i} = expression;
    expression.value = expression.value * coeffTbl.(conditions{i})(strcmp(coeffTbl.Metabolite,'Protein'));
    master_expression_perDW{i} = expression;
end
master_expression_pro_perPro = master_expression_perPro;
master_expression_pro_perDW = master_expression_perDW;

% load microarray ==> we dont analyze RNA data since it is not clean
% rnaTbl = readtable('./../input/YeastJoshua/MicroArray/matched_knn_imputed_log2_FC_to_reference_pmid_17959824.txt');% this is the log2(FC_reference)
% we use the non-log level
% rnaTbl{:,3:end} = 2.^rnaTbl{:,3:end};

% preprocess the expression table
% to facilate the future use of the expression of many samples, we
% re-organize it into a structure variable.
% the FPA matrix will be in the same order as the master_expression
% make a new master_expression for these four conditions.
% master_expression_perRNA = {};% we call this variable "master_expression"
% master_expression_perDW= {};% we call this variable "master_expression"
% geneInd = ismember(rnaTbl.YORF, model.genes); % get the index of genes in the model
% for i = 1:length(conditions)
%     expression = struct();
%     expression.genes = rnaTbl.YORF(geneInd);
%     expression.value = rnaTbl.(conditions{i})(geneInd);
%     master_expression_perRNA{i} = expression;
%     expression.value = expression.value * coeffTbl.(conditions{i})(strcmp(coeffTbl.Metabolite,'RNA'));
%     master_expression_perDW{i} = expression;
% end
% master_expression_rna = master_expression_perRNA;
% master_expression_rna_perDW = master_expression_perDW;
%% parpool
parpool(4)
%% calculate penalty (normalized level)
penalty = calculatePenalty_partialExcluded(model,master_expression_pro_perPro);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only data-derived level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns = model.rxns(~rmInd);
% back up
normalizedLevel_pro_perPro = normalizedLevel;
valid_rxns_pro_perPro = valid_rxns;

%% calculate penalty (normalized level)
penalty = calculatePenalty_partialExcluded(model,master_expression_pro_perDW);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only data-derived level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns = model.rxns(~rmInd);
% back up
normalizedLevel_pro_perDW = normalizedLevel;
valid_rxns_pro_perDW = valid_rxns;
%%
save('output/normalizedLevels_partialExcluded.mat','normalizedLevel_pro_perPro','valid_rxns_pro_perPro','normalizedLevel_pro_perDW','valid_rxns_pro_perDW');
%% calculate penalty (normalized level)
penalty = calculatePenalty_singlePartialExcluded(model,master_expression_pro_perPro);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only data-derived level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns = model.rxns(~rmInd);
% back up
normalizedLevel_pro_perPro = normalizedLevel;
valid_rxns_pro_perPro = valid_rxns;

%% calculate penalty (normalized level)
penalty = calculatePenalty_singlePartialExcluded(model,master_expression_pro_perDW);
normalizedLevel = 1 ./ penalty;
normalizedLevel(:,end) = [];
% only data-derived level will be tested
rmInd = all(normalizedLevel == 1,2);
normalizedLevel(rmInd,:) = [];
valid_rxns = model.rxns(~rmInd);
% back up
normalizedLevel_pro_perDW = normalizedLevel;
valid_rxns_pro_perDW = valid_rxns;
%%
save('output/normalizedLevels_singlePartialExcluded.mat','normalizedLevel_pro_perPro','valid_rxns_pro_perPro','normalizedLevel_pro_perDW','valid_rxns_pro_perDW');
%%
treatments = {'','_singlePartialExcluded','_partialExcluded'};
%%
for kk = 1:length(treatments)
    close all
    treatment = treatments{kk};
    load(['output/normalizedLevels',treatment,'.mat'])
    %% correlation with immediate vincinity: protein : real unit
    normalizedLevel = normalizedLevel_pro_perDW;
    valid_rxns = valid_rxns_pro_perDW ;
    % flux table
    fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
    % make the matched matrix
    fluxMat =[];
    for i = 1: length(conditions)
        fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
    end
    % normalize flux unit
    % flux is in  (moles / hr / mL cells); could be further normalized to
    % mmole/hr/gDW by chemostat info: gDCW/ml. such that it is comparable with
    % the per gDW protein content
    dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
    factor = repmat(dwTbl.gDCW_mL',size(fluxMat,1),1);
    fluxMat = fluxMat * 1000 ./ factor; %mmoles/hr/gDW
    fluxMat_raw = fluxMat;

    %Computing the correlation between a reaction expression and measured growth rate
    rho=[];
    r = [];
    p_rho=[];
    p_r=[];
    testedRxn = {};
    rxnLabel = fluxTbl.Model_Reaction_ID;
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            [rho(end+1),p_rho(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
            [r(end+1),p_r(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        end
    end
    % correct multiple testing 
    fdr_rho = mafdr(p_rho,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))

    fprintf('%d rxns give significant positive correlation by spearman\n',sum(rho(fdr_rho<0.05)>0));
    fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

    sigCorrRxns = testedRxn(fdr_r<0.05 & r >0);
    corrTbl = table(testedRxn',r',fdr_r');
    save(['output/partialGPR/',treatment,'sig_corr_rxns_raw2raw.mat'],'sigCorrRxns','corrTbl');

    figure(1)
    hold on
    histogram(rho,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
    xlim([-1,1]);
    xlabel('Correlation coefficient');
    ylabel('Number of reactions');
    histogram(rho(fdr_rho<0.05 & rho >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
    legend({'all testable reactions','significantly correlated reactions'})
    hold off
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [5.2, 4.3];
    plt.LineWidth = 2;
    plt.FontSize = 15;
    plt.XTick = -1:0.2:1;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.export(['figures/partialGPR/',treatment,'Expression_flux_correlation_rawFlux_spearman.tiff']);

    figure(2)
    hold on
    histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
    xlim([-1,1]);
    xlabel('Correlation coefficient');
    ylabel('Number of reactions');
    histogram(r(fdr_r<0.05 & r >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
    legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
    hold off
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [5.2, 4.3];
    plt.LineWidth = 2;
    plt.FontSize = 15;
    plt.XTick = -1:0.2:1;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.export(['figures/partialGPR/',treatment,'Expression_flux_correlation_rawFlux_pearson.tiff']);
    
    %% check enrichment for inrreversibe rxns
    Candidates = testedRxn(fdr_r<0.05 & r >0);
    candidatesP = fdr_r(fdr_r<0.05 & r >0);
    [~,A] = sort(candidatesP);
    Candidates = Candidates(A);
    
    hitRxns = Candidates;
    allRxnSet_tested = testedRxn;
    allRxnSet_GPR = model.rxns(cellfun(@(x) ~isempty(x),model.grRules));

    N_irr_hitRxns = sum(~(model.lb(ismember(model.rxns,hitRxns)) < -1e-10 & model.ub(ismember(model.rxns,hitRxns)) > 1e-10));
    N_irr_allRxnSet_tested = sum(~(model.lb(ismember(model.rxns,allRxnSet_tested)) < -1e-10 & model.ub(ismember(model.rxns,allRxnSet_tested)) > 1e-10));
    N_irr_allRxnSet_GPR = sum(~(model.lb(ismember(model.rxns,allRxnSet_GPR)) < -1e-10 & model.ub(ismember(model.rxns,allRxnSet_GPR)) > 1e-10));

    N_irr_hitRxns
    length(allRxnSet_tested)
    N_irr_allRxnSet_tested
    length(hitRxns)
    p1 = 1-hygecdf(N_irr_hitRxns-1,length(allRxnSet_tested),N_irr_allRxnSet_tested,length(hitRxns))% enrichment in all tested rxns
    p2 = 1-hygecdf(N_irr_hitRxns-1,length(allRxnSet_GPR),N_irr_allRxnSet_GPR,length(hitRxns))% enrichment in all rxns (with GPR)

    %% correlation with immediate vincinity: protein : normalized unit free
    normalizedLevel = normalizedLevel_pro_perPro;
    valid_rxns = valid_rxns_pro_perPro ;
    % flux table 
    fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
    % make the matched matrix
    fluxMat = [];
    for i = 1: length(conditions)
        fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
    end
    % normalize flux unit
    % when normalize the flux by the flux of biomass production (growth rate),
    % the unit of growth rate needs to be taken care of. In chemostat setting,
    % steady state was defined as stable OD (see SIMMER paper), which means
    % steady cell density (number, aka, volume). Therefore, the dilution rate
    % is a measure of per cell flux. So, we should normalize the internal flux
    % under /ml cell metric

    % flux is in  (moles / hr / mL cells); no conversion is needed. 
    % in fact, correlation got worse if we normalzie the flux to / gDW first!

    GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
    fluxMat = fluxMat ./ repmat(GRrate.DR_Actual',size(fluxMat,1),1);
    fluxMat_normalized = fluxMat;
    %Computing the correlation between a reaction expression and measured growth rate
    rho=[];
    p_rho=[];
    r=[];
    p_r=[];
    testedRxn = {};
    rxnLabel = fluxTbl.Model_Reaction_ID;
    for j = 1:length(rxnLabel)
        fluxMeasure = fluxMat(j,:);
        if any(strcmp(valid_rxns,rxnLabel{j}))
            testedRxn(end+1) = rxnLabel(j);
            [rho(end+1),p_rho(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
            [r(end+1),p_r(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
            %Correcting for multiple hypothesis using FDR and significance
            %level of?
        end
    end
    fdr_rho = mafdr(p_rho,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
    fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))

    fprintf('%d rxns give significant positive correlation by spearman\n',sum(rho(fdr_rho<0.05)>0));
    fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

    sigCorrRxns = testedRxn(fdr_r<0.05 & r >0);
    corrTbl = table(testedRxn',r',fdr_r');
    save(['output/partialGPR/',treatment,'sig_corr_rxns_rel2rel.mat'],'sigCorrRxns','corrTbl');

    figure(4)
    hold on
    histogram(rho,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
    xlim([-1,1]);
    xlabel('Correlation coefficient');
    ylabel('Number of reactions');
    histogram(rho(fdr_rho<0.05 & rho >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
    legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
    hold off
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [5.2, 4.3];
    plt.LineWidth = 2;
    plt.FontSize = 15;
    plt.XTick = -1:0.2:1;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.export(['figures/partialGPR/',treatment,'Expression_flux_correlation_relFlux_spearman.tiff']);

    figure(5)
    hold on
    histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
    xlim([-1,1]);
    xlabel('Correlation coefficient');
    ylabel('Number of reactions');
    histogram(r(fdr_r<0.05 & r >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
    legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
    hold off
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [5.2, 4.3];
    plt.LineWidth = 2;
    plt.FontSize = 15;
    plt.XTick = -1:0.2:1;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.export(['figures/partialGPR/',treatment,'Expression_flux_correlation_relFlux_pearson.tiff']);
    %% individual plot
    Candidates = testedRxn(fdr_r<0.05 & r >0);
    candidatesP = fdr_r(fdr_r<0.05 & r >0);
    [~,A] = sort(candidatesP);
    Candidates = Candidates(A);
    %% check enrichment for inrreversibe rxns
    hitRxns = Candidates;
    allRxnSet_tested = testedRxn;
    allRxnSet_GPR = model.rxns(cellfun(@(x) ~isempty(x),model.grRules));

    N_irr_hitRxns = sum(~(model.lb(ismember(model.rxns,hitRxns)) < -1e-10 & model.ub(ismember(model.rxns,hitRxns)) > 1e-10));
    N_irr_allRxnSet_tested = sum(~(model.lb(ismember(model.rxns,allRxnSet_tested)) < -1e-10 & model.ub(ismember(model.rxns,allRxnSet_tested)) > 1e-10));
    N_irr_allRxnSet_GPR = sum(~(model.lb(ismember(model.rxns,allRxnSet_GPR)) < -1e-10 & model.ub(ismember(model.rxns,allRxnSet_GPR)) > 1e-10));

    N_irr_hitRxns
    length(allRxnSet_tested)
    N_irr_allRxnSet_tested
    length(hitRxns)
    p1 = 1-hygecdf(N_irr_hitRxns-1,length(allRxnSet_tested),N_irr_allRxnSet_tested,length(hitRxns))
    p2 = 1-hygecdf(N_irr_hitRxns-1,length(allRxnSet_GPR),N_irr_allRxnSet_GPR,length(hitRxns))
    p3 = 1-hygecdf(N_irr_allRxnSet_tested-1,length(allRxnSet_GPR),N_irr_allRxnSet_GPR,length(allRxnSet_tested)) % the tested rxns are already enriched 

    %% enrichment of direct connections
    byProducts = {'carbon dioxide';'AMP';'NADP(+)';'NADPH';'diphosphate';'oxygen';'NADH';'NAD';'phosphate';'ADP';'coenzyme A';'ATP';'H2O';'H+';'GTP';'GDP';'(R)-carnitine';'FAD';'FADH2'};
    % add compartment label to byproducts
    byProducts = model.mets(ismember(cellfun(@(x) regexprep(x,'\s\[.*\]$',''),model.metNames, 'UniformOutput',false),byProducts));
    model_bak = model;
    model.S(ismember(model.mets,byProducts),:) = [];
    model.mets(ismember(model.mets,byProducts)) = [];

    cMat = zeros(length(model.rxns),length(model.rxns));% connection matrix
    for i = 1:length(model.rxns)
        % out degree
        if model.lb(i) < -1e-10
            metsID = model.S(:,i)<0;
            nextRxnID = any(model.S(metsID,:),1);
            nextRxnID(i) = false;
            rxnset = model.rxns(nextRxnID);
            r1 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)>0,1)');...
                rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)<0,1)')];
        else
            r1 = {};
        end
        if model.ub(i) > 1e-10
            metsID = model.S(:,i)>0;
            nextRxnID = any(model.S(metsID,:),1);
            nextRxnID(i) = false;
            rxnset = model.rxns(nextRxnID);
            r2 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)>0,1)');...
                rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)<0,1)')];
        else
            r2 = {};
        end

        % in degree
        if model.lb(i) < -1e-10
            metsID = model.S(:,i)>0;
            nextRxnID = any(model.S(metsID,:),1);
            nextRxnID(i) = false;
            rxnset = model.rxns(nextRxnID);
            r3 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)<0,1)');...
                rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)>0,1)')];
        else
            r3 = {};
        end
        if model.ub(i) > 1e-10
            metsID = model.S(:,i)<0;
            nextRxnID = any(model.S(metsID,:),1);
            nextRxnID(i) = false;
            rxnset = model.rxns(nextRxnID);
            r4 = [rxnset(model.lb(nextRxnID) < -1e-10 & any(model.S(metsID,nextRxnID)<0,1)');...
                rxnset(model.ub(nextRxnID) > 1e-10 & any(model.S(metsID,nextRxnID)>0,1)')];
        else
            r4 = {};
        end
        r_c = unique([r1;r2;r3;r4]);
        cMat(i,ismember(model.rxns,r_c)) = 1;
    end
    model = model_bak;       

    % analysis 
    hitRxns = Candidates;
    allRxnSet_tested = testedRxn;
    allRxnSet_GPR = model.rxns(cellfun(@(x) ~isempty(x),model.grRules));

    c_hit = sum(sum(cMat(ismember(model.rxns,hitRxns),ismember(model.rxns,hitRxns))));
    % permutation test - against all tested rxns
    c_perm = zeros(100000,1);
    rng(1126)
    for i = 1: 100000
        sample = allRxnSet_tested(randperm(length(allRxnSet_tested),length(hitRxns)));
        c_perm(i) = sum(sum(cMat(ismember(model.rxns,sample),ismember(model.rxns,sample))));
    end

    figure;
    histogram(c_perm,'Normalization','pdf','BinEdges',0:2:140)
    xline(c_hit,'r')
    xlabel('Total number of connections','FontSize', 12)
    ylabel('Frequency','FontSize', 12)
    txt = ['p = ', num2str((1+ sum(c_perm>=c_hit)) / (1+length(c_perm)))];
    text(c_hit+5,0.025,txt,'HorizontalAlignment','left','FontSize', 15)
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [5.2, 4.3];
    plt.LineWidth = 2;
    plt.FontSize = 15;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.export(['figures/partialGPR/',treatment,'Expression_flux_correlation_connectivityEnrichment.tiff']);
end
% %% permutation of labels to show if these 35% enrichment is significant
% % permute y
% percSig_sample = [];
% percSig_posCorr_sample = [];
% percSig_posCorr_sample2 = [];
% rng(1126)
% nperm = 10000;
% for i = 1: nperm
%     %Computing the correlation between a reaction expression and measured growth rate
%     rho=[];
%     p=[];
%     testedRxn = {};
%     rxnLabel = fluxTbl.Model_Reaction_ID;
%     fluxMat_perm = fluxMat(:,randperm(size(fluxMat,2)));
%     for j = 1:length(rxnLabel)
%         fluxMeasure = fluxMat_perm(j,:);
%         if any(strcmp(valid_rxns,rxnLabel{j}))
%             testedRxn(end+1) = rxnLabel(j);
%             [rho(end+1),p(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
%             %Correcting for multiple hypothesis using FDR and significance level of
%         end
%     end
%     percSig_sample = [percSig_sample,sum(p<0.05)/length(testedRxn)];
%     percSig_posCorr_sample = [percSig_posCorr_sample,sum(rho(p<0.05)>0)/length(rho(p<0.05))];
%     percSig_posCorr_sample2 = [percSig_posCorr_sample2,sum(rho(p<0.05)>0)/length(testedRxn)];
% end
% percSig_posCorr_sample(isnan(percSig_posCorr_sample)) = 0;
% figure(1)
% histogram(percSig_sample)
% sum(percSig_sample >= 0.3554)/nperm%perm 60,000 p= 0.0028
% figure(2)
% histogram(percSig_posCorr_sample)
% sum(percSig_posCorr_sample >= 0.7288)/nperm
% figure(3)
% histogram(percSig_posCorr_sample2)
% sum(percSig_posCorr_sample2 >= 0.2590)/nperm%perm 60,000 p= 0.0025
% %% permute x ==> same result
% percSig_sample = [];
% percSig_posCorr_sample = [];
% percSig_posCorr_sample2 = [];
% rng(1126)
% nperm = 1000;
% for i = 1: nperm
%     %Computing the correlation between a reaction expression and measured growth rate
%     rho=[];
%     p=[];
%     testedRxn = {};
%     rxnLabel = fluxTbl.Model_Reaction_ID;
%     normalizedLevel_perm = normalizedLevel(:,randperm(size(normalizedLevel,2)));
%     for j = 1:length(rxnLabel)
%         fluxMeasure = fluxMat(j,:);
%         if any(strcmp(valid_rxns,rxnLabel{j}))
%             testedRxn(end+1) = rxnLabel(j);
%             [rho(end+1),p(end+1)] = corr(normalizedLevel_perm(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
%             %Correcting for multiple hypothesis using FDR and significance level of
%         end
%     end
%     percSig_sample = [percSig_sample,sum(p<0.05)/length(testedRxn)];
%     percSig_posCorr_sample = [percSig_posCorr_sample,sum(rho(p<0.05)>0)/length(rho(p<0.05))];
%     percSig_posCorr_sample2 = [percSig_posCorr_sample2,sum(rho(p<0.05)>0)/length(testedRxn)];
% end
% percSig_posCorr_sample(isnan(percSig_posCorr_sample)) = 0;
% figure(1)
% histogram(percSig_sample)
% sum(percSig_sample >= 0.3554)/nperm
% figure(2)
% histogram(percSig_posCorr_sample)
% sum(percSig_posCorr_sample >= 0.7288)/nperm
% figure(3)
% histogram(percSig_posCorr_sample2)
% sum(percSig_posCorr_sample2 >= 0.2590)/nperm
% %% correlation with any reaction: protien
% normalizedLevel = normalizedLevel_pro;
% valid_rxns = valid_rxns_pro ;
% % flux table
% fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% % make the matched matrix
% for i = 1: length(conditions)
%     fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
% end
% %Computing the correlation between a reaction expression and measured growth rate
% sigLevels={};
% sigRho={};
% sigRxns={};
% sigFluxMeas={};
% rxnLabel = fluxTbl.Model_Reaction_ID;
% for j = 1:length(rxnLabel)
%     fluxMeasure = fluxMat(j,:);
%     rho = [];
%     p = [];
%     for i=1:size(normalizedLevel,1)
%         [rho(i),p(i)] = corr(normalizedLevel(i,:)',abs(fluxMeasure)','type','Spearman');
%     end
%     %Correcting for multiple hypothesis using FDR and significance level of
%     alpha = 0.05;
%     [pID,pN] = FDR(p,alpha);
%     if ~isempty(pID)
%         fprintf('sig corr found for reaction %s (%s)\n',rxnLabel{j},model.rxnNames{strcmp(model.rxns,rxnLabel{j})});
% 
%         %Identifying the set of significantly correlated reactions
%         [p_sort,idx_sort] = sort(p);
%         sig_p = find(p_sort<=pID);
%         idx = idx_sort(1:length(sig_p));
% 
%         %Update all variables to include only the set found above
%         sigLevels{end+1} = normalizedLevel(idx,:);
%         sigRho{end+1} = rho(idx);
%         sigRxns{end+1} = valid_rxns(idx);
%         sigFluxMeas{end+1} = rxnLabel{j};
%     end
% end
% sigLevels_pro = sigLevels;
% sigRho_pro = sigRho;
% sigRxns_pro = sigRxns;
% sigFluxMeas_pro = sigFluxMeas;
% %% plot example 
% fluxQry = 'r_0959';%r_0816 r_0959 r_0152 r_0727 r_0889
% expTarget = 'r_0959';
% fluxMeasure = fluxMat(strcmp(rxnLabel,fluxQry),:);
% corr(normalizedLevel(strcmp(valid_rxns,expTarget),:)',abs(fluxMeasure)','type','Spearman')
% fit = fitlm(normalizedLevel(strcmp(valid_rxns,expTarget),:)',abs(fluxMeasure)');
% plot(fit)
% %% add FVA errors 
% for i = 1: length(conditions)
%     errMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) - fluxTbl.([conditions{i},'_FVA_MIN']))/2;
% end
% errMeasure = errMat(strcmp(rxnLabel,fluxQry),:);
% hold on 
% errorbar(normalizedLevel(strcmp(valid_rxns,expTarget),:)',abs(fluxMeasure)',errMeasure,'o');
% hold off
%% ===============================
%% ===============================
%% ===============================
% %% compare the agreement of protein and mRNA
% target = intersect(intersect(rnaTbl.YORF,proTbl.Gene),model.genes);
% vector_pro = [];
% vector_rna = [];
% for j = 1:length(target)
%     tmp_pro = [];
%     tmp_rna = [];
%     for i = 1:length(conditions)
%         [A B]= ismember(target{j},proTbl.Gene);
%         tmp_pro = [tmp_pro;proTbl.(conditions{i})(B(A))];
%         [A B]= ismember(target{j},rnaTbl.YORF);
%         tmp_rna = [tmp_rna;rnaTbl.(conditions{i})(B(A))];
%     end
%     vector_pro = [vector_pro;tmp_pro./tmp_pro(1)];
%     vector_rna = [vector_rna;tmp_rna./tmp_pro(1)];
% end
% figure
% fit = fitlm(log2(vector_rna),log2(vector_pro));
% plot(fit)
