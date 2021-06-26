addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/

%% 1. load the model and prepare the model
addpath('./../scripts/')
model = loadYeatModel();
% the following nutrients need to be set manually
% to start with basic FPA, we allow unlimited exchange
% phosphate exchange
model.lb(strcmp(model.rxns,'r_2005')) = -1000;
% glucose exchange
model.lb(strcmp(model.rxns,'r_1714')) = -1000;
% ammonium exchange 
model.lb(strcmp(model.rxns,'r_1654')) = -1000;
% uracil 
model.lb(strcmp(model.rxns,'r_2090')) = -1000;
% leucine
model.lb(strcmp(model.rxns,'r_1899')) = -1000;
% maintanence 
model = changeRxnBounds(model,'r_4046',0,'l'); % maintance 
model = changeRxnBounds(model,'r_4046',1000,'u'); % maintance 

% decouple energy from growth
% model.S(ismember(model.mets,{'s_0434[c]','s_0803[c]','s_0394[c]','s_0794[c]','s_1322[c]'}), strcmp('r_4041',model.rxns)) = 0; 

% allow all exchange
% model.lb(findExcRxns(model)) = -1000;
%% 3. load the expression files, distance matrix, and other optional inputs
% load proteomics
% proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
% % we use the non-log level
% proTbl{:,2:end} = 2.^proTbl{:,2:end};
% % preprocess the expression table
% % to facilate the future use of the expression of many samples, we
% % re-organize it into a structure variable.
% % the FPA matrix will be in the same order as the master_expression
% conditions = proTbl.Properties.VariableNames(2:end);
% % make a new master_expression for these four conditions.
% master_expression = {};% we call this variable "master_expression"
% geneInd = ismember(proTbl.Gene, model.genes); % get the index of genes in the model
% for i = 1:length(conditions)
%     expression = struct();
%     expression.genes = proTbl.Gene(geneInd);
%     expression.value = proTbl.(conditions{i})(geneInd);
%     master_expression{i} = expression;
% end
% master_expression_pro = master_expression;
% 
% % load microarray
% rnaTbl = readtable('./../input/YeastJoshua/MicroArray/matched_knn_imputed_log2_FC_to_reference_pmid_17959824.txt');% this is the log2(FC_reference)
% % we use the non-log level
% rnaTbl{:,3:end} = 2.^rnaTbl{:,3:end};
% % preprocess the expression table
% % to facilate the future use of the expression of many samples, we
% % re-organize it into a structure variable.
% % the FPA matrix will be in the same order as the master_expression
% % make a new master_expression for these four conditions.
% master_expression = {};% we call this variable "master_expression"
% geneInd = ismember(rnaTbl.YORF, model.genes); % get the index of genes in the model
% for i = 1:length(conditions)
%     expression = struct();
%     expression.genes = rnaTbl.YORF(geneInd);
%     expression.value = rnaTbl.(conditions{i})(geneInd);
%     master_expression{i} = expression;
% end
% master_expression_rna = master_expression;

% load the distance matrix
% users can uncomment the following codes to load from distance calculator
% output; here we load directly from saved matlab variable because of file
% size restriction of GitHub
distance_raw = readtable('./../input/YeastJoshua/distanceMatrix_weighted.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github
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
% load the special penalties 
% In general, we recommend tp set penalty for all Exchange, Demand,
% and Sink reactions to 0 to not penaltize the external reactions. Users 
% may need to interactively tune their special penalties for best flux
% distribution in the FPA calculation
extRxns = model.rxns(findExcRxns(model));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));

% we dont recomand any specific special distance for generic model; In the
% dual model, the special distance was used to discourage using of side
% metabolites. Since side/storage metabolites are not applicable for
% generic model, we don't use any special distance. 
% manualDist = {};

proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

% make the blocklist to represent the condition-specific nutrient 
p_lim = {'r_1899_r','r_2090_r'};
c_lim = {'r_1899_r','r_2090_r'};
n_lim = {'r_1899_r','r_2090_r'};
l_lim = {'r_2090_r'};
u_lim = {'r_1899_r'};
blocklist = [repmat({p_lim},1,5),repmat({c_lim},1,5),repmat({n_lim},1,5),repmat({l_lim},1,5),repmat({u_lim},1,5),{{}}];
%% flux table 
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;
% normalize flux unit
% when normalize the flux by the flux of biomass production (growth rate),
% the unit of growth rate needs to be taken care of. In chemostat setting,
% steady state was defined as stable OD (see SIMMER paper), which means
% steady cell density (number, aka, volume). Therefore, the dilution rate
% is a measure of per cell flux. So, we should normalize the internal flux
% under /ml cell metric

% flux is in  (moles / hr / mL cells); no conversion is needed. 
% in fact, correlation got worse if we normalzie the flux to / gDW first!
fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
%  dilutionFactor = repmat([0.05 0.1 0.16 0.22 0.30],size(fluxMat,1),5);
%  fluxMat_double_normalized = fluxMat_normalized ./ dilutionFactor;

% make the raw flux matrix (in per gDW unit)
% flux is in  (moles / hr / mL cells); could be further normalized to
% mmole/hr/gDW by chemostat info: gDCW/ml. such that it is comparable with
% the per gDW protein content
dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
fluxMat_raw = fluxMat;
factor = repmat(dwTbl.gDCW_mL',size(fluxMat_raw,1),1);
fluxMat_raw = fluxMat_raw * 1000 ./ factor; %mmoles/hr/gDW

%% find internal reactions
% EXrelRxns = model.rxns(any(model.S(cellfun(@(x) ~isempty(regexp(x,'\[e\]$','once')),model.mets),:),1));
% transporters = intersect(rxnLabel,EXrelRxns);
% printRxnFormula_XL(model,transporters);
% internalRxnInd = ~ismember(rxnLabel,transporters);
%% load merged levels
load('output/normalizedLevels_partialExcluded.mat');
%% precalc penalty 
% penalty_merged = ones(length(model.rxns),size(fluxMat,2)+1);
% [A B] = ismember(model.rxns,valid_rxns_merged);
% penalty_merged(A,1:(end-1)) = ones(size(normalizedLevel_merged(B(A),:),1),size(normalizedLevel_merged(B(A),:),2)) ./ normalizedLevel_merged(B(A),:);

penalty_pro = ones(length(model.rxns),size(fluxMat,2)+1);
[A B] = ismember(model.rxns,valid_rxns_pro_perPro);
penalty_pro(A,1:(end-1)) = ones(size(normalizedLevel_pro_perPro(B(A),:),1),size(normalizedLevel_pro_perPro(B(A),:),2)) ./ normalizedLevel_pro_perPro(B(A),:);

penalty_pro_raw = ones(length(model.rxns),size(fluxMat,2)+1);
[A B] = ismember(model.rxns,valid_rxns_pro_perDW);
penalty_pro_raw(A,1:(end-1)) = ones(size(normalizedLevel_pro_perDW(B(A),:),1),size(normalizedLevel_pro_perDW(B(A),:),2)) ./ normalizedLevel_pro_perDW(B(A),:);

% apply additional penalty to the exchange of limiting nutrients
for i = 1:size(manualPenalty,1)
    penalty_pro(strcmp(model.rxns,manualPenalty{i,1}),:) = manualPenalty{i,2};
end
penalty_pro(strcmp(model.rxns,'r_2005'),1:5) = 10;
penalty_pro(strcmp(model.rxns,'r_1714'),6:10) = 10;
penalty_pro(strcmp(model.rxns,'r_1654'),11:15) = 10;

for i = 1:size(manualPenalty,1)
    penalty_pro_raw(strcmp(model.rxns,manualPenalty{i,1}),:) = manualPenalty{i,2};
end
penalty_pro_raw(strcmp(model.rxns,'r_2005'),1:5) = 10;
penalty_pro_raw(strcmp(model.rxns,'r_1714'),6:10) = 10;
penalty_pro_raw(strcmp(model.rxns,'r_1654'),11:15) = 10;
%% evaluate different FPA
setups = {'oriDist_oriDecay','wtdDist_oriDecay'};
for zz = 1:2
    %% correlation for 1-d titration
    load(['output/Titration_relativeExp_',setups{zz},'.mat'])
    dorders = n;
    rMat = zeros(length(targetRxns),length(dorders));
    N_sigCorr = zeros(length(dorders),1);
    N_sigCorr_highCV = zeros(length(dorders),1);
    perc_sigCorr_in_highCV = zeros(length(dorders),1);
    for nn = 1: length(dorders)
        %%
        % nn = 10;
        dorders(nn);
        FP = FP_collection{nn};
        relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
        for i = 1:size(FP,1)
            for j = 1:(size(FP,2)-1)
                if  mean(fluxMat(i,:)) > 0 
                    if ~isnan(FP{i,j}(1))
                        relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                    elseif ~isnan(FP{i,j}(2))
                        relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                    else
                        error('check');
                    end
                else
                    if ~isnan(FP{i,j}(2))
                        relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                    elseif ~isnan(FP{i,j}(1))
                        relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                    else
                        error('check');
                    end
                end
            end
        end
        rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
        relFP = relFP(~rmInd,:);
        valid_rxns = targetRxns(~rmInd);
        %Computing the correlation between a reaction expression and measured growth rate
        r=[];
        p_r=[];
        deltaminmax = [];
        testedRxn = {};
        for j = 1:length(rxnLabel)
            fluxMeasure = fluxMat_normalized(j,:);
            if any(strcmp(valid_rxns,rxnLabel{j}))
                testedRxn(end+1) = rxnLabel(j);
                [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
                %Correcting for multiple hypothesis using FDR and significance level of
                deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
            end
        end

        fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
        fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));

        if (dorders(nn) == 1.5)
            figure;
            hold on
            histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
            xlim([-1,1]);
            xlabel('Correlation coefficient');
            ylabel('Number of reactions');
            histogram(r(fdr_r<0.05 & r >0 & deltaminmax > 0.2),'FaceColor','#D95319','BinEdges',-1:0.2:1)
            legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
            hold off
            plt = Plot(); % create a Plot object and grab the current figure
            plt.BoxDim = [5.2, 4.3];
            plt.LineWidth = 2;
            plt.FontSize = 15;
            plt.XTick = -1:0.2:1;
            plt.LegendLoc = 'NorthWest';
            plt.FontName = 'Arial';
            plt.export(['figures/FPA_flux_correlation_',setups{zz},'_pearson.tiff']);
        end

        [A B] = ismember(testedRxn,targetRxns);
        rMat(B(A),nn) = r;
        N_sigCorr(nn) = sum(r(fdr_r<0.05)>0);
        N_sigCorr_highCV(nn) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0);
        perc_sigCorr_in_highCV(nn) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0) / sum(deltaminmax > 0.2);
    end
    %% analyze some overall metric ==> will be plotted together with others
    % figure(1)
    % plot(dorders,N_sigCorr,'.-')
    % xlabel('distance order');
    % ylabel('number of sig corr rxns');
    % figure(2)
    % plot(dorders,N_sigCorr_highCV,'.-')
    % xlabel('distance order');
    % ylabel('number of sig-corr & rFP-CV > 0.1 rxns ');
    % figure(3)
    % plot(dorders,perc_sigCorr_in_highCV,'.-')
    % xlabel('distance order');
    % ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');

    %% save data to plot four setups together: 1-d tritrations
    if zz==1
        %decayBound = exp(log(1./minInfo) ./ dorders)-1;
        %dorders_normal = decayBound;
        dorders_normal = dorders;% dorders/max(dorders);
        N_sigCorr_normal = N_sigCorr;
        N_sigCorr_highCV_normal = N_sigCorr_highCV;
        perc_sigCorr_in_highCV_normal = perc_sigCorr_in_highCV;
    elseif zz ==2
        %decayBound = exp(log(1./minInfo) ./ dorders)-1;
        %dorders_wtdDist = decayBound;
        dorders_wtdDist = dorders;% dorders/max(dorders);
        N_sigCorr_wtdDist = N_sigCorr;
        N_sigCorr_highCV_wtdDist = N_sigCorr_highCV;
        perc_sigCorr_in_highCV_wtdDist = perc_sigCorr_in_highCV;
    end
end

%% the 2-d titrition data
close all
setups = {'oriDist_expDecay','wtdDist_expDecay'};
for jj = 1:2
    %% correlation for 1-d titration
    load(['output/Titration_relativeExp_',setups{jj},'.mat'])
    %% correlation for 2-d titration
    N_sigCorr = zeros(length(n2),length(base));
    N_sigCorr_highCV = zeros(length(n2),length(base));
    perc_sigCorr_in_highCV = zeros(length(n2),length(base));
    for zz = 1:length(base)
        dorders = n2;
        rMat = zeros(length(targetRxns),length(dorders));
        for nn = 1: length(dorders)
            %%
            % nn = 10;
            %dorders(nn)
            FP = FP_collection_2{zz}{nn};
            relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
            for i = 1:size(FP,1)
                for j = 1:(size(FP,2)-1)
                    if  mean(fluxMat(i,:)) > 0 
                        if ~isnan(FP{i,j}(1))
                            relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                        elseif ~isnan(FP{i,j}(2))
                            relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                        else
                            fprintf('base = %f, n = %d, i = %d, j = %d, is a failed FPA\n',...
                                base(zz), dorders(nn),i,j);
                            warning('check');
                        end
                    else
                        if ~isnan(FP{i,j}(2))
                            relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
                        elseif ~isnan(FP{i,j}(1))
                            relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
                        else
                            error('check');
                        end
                    end
                end
            end
            rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
            relFP = relFP(~rmInd,:);
            valid_rxns = targetRxns(~rmInd);
            %Computing the correlation between a reaction expression and measured growth rate
            r=[];
            p_r=[];
            deltaminmax = [];
            testedRxn = {};
            for j = 1:length(rxnLabel)
                fluxMeasure = fluxMat_normalized(j,:);
                if any(strcmp(valid_rxns,rxnLabel{j}))
                    testedRxn(end+1) = rxnLabel(j);
                    [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
                    %Correcting for multiple hypothesis using FDR and significance level of
                    deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
                end
            end
            
            fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
            fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05)>0));
      
            if (dorders(nn) == 6 && base(zz) ==2)
                figure;
                hold on
                histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
                xlim([-1,1]);
                xlabel('Correlation coefficient');
                ylabel('Number of reactions');
                histogram(r(fdr_r<0.05 & r >0 & deltaminmax > 0.2),'FaceColor','#D95319','BinEdges',-1:0.2:1)
                legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
                hold off
                plt = Plot(); % create a Plot object and grab the current figure
                plt.BoxDim = [5.2, 4.3];
                plt.LineWidth = 2;
                plt.FontSize = 15;
                plt.XTick = -1:0.2:1;
                plt.LegendLoc = 'NorthWest';
                plt.FontName = 'Arial';
                plt.export(['figures/FPA_flux_correlation_',setups{jj},'_pearson.tiff']);
            end

            [A B] = ismember(testedRxn,targetRxns);
            rMat(B(A),nn) = r;
            N_sigCorr(nn,zz) = sum(r(fdr_r<0.05)>0);
            N_sigCorr_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0);
            perc_sigCorr_in_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0) / sum(deltaminmax > 0.2);
        end
    end
    %% analyze some overall metric
%     figure(4)
%     hold on 
%     for zz = 1:length(base)
%         plot(dorders,N_sigCorr(:,zz),'.-')
%     end
%     hold off
%     xlabel('distance order');
%     ylabel('number of sig corr rxns');
%     legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))
    colorList = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k'};
    figure;
    hold on 
    for zz = 1:length(base)
        plot(dorders,N_sigCorr_highCV(:,zz),'.-','Color',colorList{zz})
    end
    hold off
    xlabel('Distance order');
    ylabel('Number of significantly correlated reactions ');
    legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [5.2, 4.3];
    plt.LineWidth = 2;
    plt.FontSize = 15;
    plt.LegendLoc = 'NorthEast';
    plt.FontName = 'Arial';
    plt.export(['figures/FPAtitrationPlots/Nsig_vs_distanceOrder',setups{jj},'.tiff']);
%     figure(6)
%     hold on 
%     for zz = 1:length(base)
%         plot(dorders,perc_sigCorr_in_highCV(:,zz),'.-')
%     end
%     hold off
%     xlabel('distance order');
%     ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');
%     legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))

    figure;
    hold on 
    for zz = 1:length(base)
        plot(dorders(1:13),N_sigCorr_highCV(1:13,zz),'.-','Color',colorList{zz})
    end
    hold off
    xlabel('Distance order');
    ylabel('Number of significantly correlated reactions');
    legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [5.2, 4.3];
    plt.LineWidth = 2;
    plt.FontSize = 15;
    if jj ==1
        plt.LegendLoc = 'NorthEast';
    elseif jj ==2
        plt.LegendLoc = 'South';
    end
    plt.FontName = 'Arial';
    plt.export(['figures/FPAtitrationPlots/Nsig_vs_distanceOrder',setups{jj},'_zoomIn.tiff']);
    
     if jj==1
        %% exp decay only -2d
        %dorders_expDecay100 = log(1./minInfo-1) ./ log(100) + dorders;
        dorders_expDecay100 = dorders;% 1-dorders/max(dorders);
        N_sigCorr_expDecay100 = N_sigCorr(:,7);% base =100
        N_sigCorr_highCV_expDecay100 = N_sigCorr_highCV(:,7);
        perc_sigCorr_in_highCV_expDecay100 = perc_sigCorr_in_highCV(:,7);

        %dorders_expDecay2 = log(1./minInfo-1) ./ log(2) + dorders;
        dorders_expDecay2 = dorders;% 1-dorders/max(dorders);
        N_sigCorr_expDecay2 = N_sigCorr(:,2);% base =2
        N_sigCorr_highCV_expDecay2 = N_sigCorr_highCV(:,2);
        perc_sigCorr_in_highCV_expDecay2 = perc_sigCorr_in_highCV(:,2);
     elseif jj ==2
    %% exp decay + wtd dist -2d
        %dorders_wtdDist_expDecay100 = log(1./minInfo-1) ./ log(100) + dorders;
        dorders_wtdDist_expDecay100 = dorders;% 1-dorders/max(dorders);
        N_sigCorr_wtdDist_expDecay100 = N_sigCorr(:,7);% base =100
        N_sigCorr_highCV_wtdDist_expDecay100 = N_sigCorr_highCV(:,7);
        perc_sigCorr_in_highCV_wtdDist_expDecay100 = perc_sigCorr_in_highCV(:,7);

        %dorders_wtdDist_expDecay2 = log(1./minInfo-1) ./ log(2) + dorders;
        dorders_wtdDist_expDecay2 = dorders;% 1-dorders/max(dorders);
        N_sigCorr_wtdDist_expDecay2 = N_sigCorr(:,2);% base =2
        N_sigCorr_highCV_wtdDist_expDecay2 = N_sigCorr_highCV(:,2);
        perc_sigCorr_in_highCV_wtdDist_expDecay2 = perc_sigCorr_in_highCV(:,2);
     end
end

%% plot 
% figure(4)
% hold on
% plot(dorders_normal,N_sigCorr_normal,'.--','LineWidth',2)
% plot(dorders_wtdDist,N_sigCorr_wtdDist,'.-','LineWidth',2)
% plot(dorders_expDecay100,N_sigCorr_expDecay100,'.--','LineWidth',2)
% plot(dorders_expDecay2,N_sigCorr_expDecay2,'.--','LineWidth',2)
% plot(dorders_wtdDist_expDecay100,N_sigCorr_wtdDist_expDecay100,'.-','LineWidth',2)
% plot(dorders_wtdDist_expDecay2,N_sigCorr_wtdDist_expDecay2,'.-','LineWidth',2)
% xlabel('most global --> most local (scale not comparable!)');
% ylabel('number of sig corr rxns');
% legend({'oriDist + oriDecay','wtdDist + oriDecay',...
%     'oriDist + wtdDecay (100)','oriDist + wtdDecay (2)'...
%     'wtdDist + wtdDecay (100)','wtdDist + wtdDecay (2)'},'Location','north','FontSize',10,'FontWeight','bold');
baseline = 46;

figure;
hold on
plot(dorders_wtdDist_expDecay2(1:13),N_sigCorr_highCV_wtdDist_expDecay2(1:13) ,'o-','LineWidth',2,'Color','#D95319')
plot(dorders_wtdDist_expDecay100(1:13),N_sigCorr_highCV_wtdDist_expDecay100(1:13) ,'o-','LineWidth',2,'Color','#EDB120')
plot(dorders_wtdDist(1:16),N_sigCorr_highCV_wtdDist(1:16) ,'o-','LineWidth',2,'Color','#0072BD')
plot(dorders_normal(1:16),N_sigCorr_highCV_normal(1:16) ,'o--k','LineWidth',2)
% plot(dorders_expDecay100(1:13),N_sigCorr_highCV_expDecay100(1:13) ,'.--','LineWidth',2)
% plot(dorders_expDecay2(1:13),N_sigCorr_highCV_expDecay2(1:13) ,'.--','LineWidth',2)
yline(baseline,'--','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Distance order/boundary');
ylabel('Number of significantly correlated reactions ');
ylim([20 80])
legend({'improved FPA (base = 2)','improved FPA (base = 100)','original FPA + weighted distance', 'original FPA','expression only'},'Location','southeast','FontSize',12);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.7, 4.7];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/FPAbenchmarking.tiff']);

figure;
hold on
plot(dorders_wtdDist_expDecay100(1:14),N_sigCorr_highCV_wtdDist_expDecay100(1:14) ,'o-','LineWidth',2,'Color','#D95319')
% plot(dorders_expDecay100(1:13),N_sigCorr_highCV_expDecay100(1:13) ,'.--','LineWidth',2)
% plot(dorders_expDecay2(1:13),N_sigCorr_highCV_expDecay2(1:13) ,'.--','LineWidth',2)
yline(baseline,'--','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Integration boundary (+/- # of reactions)');
ylabel('Number of significantly correlated reactions ');
ylim([20 80])
legend({'network integration','expression only'},'Location','northeast','FontSize',15);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.7, 4.7];
plt.LineWidth = 2;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/FPA_localOverallBenefit.tiff']);
% figure(6)
% hold on
% plot(dorders_normal,perc_sigCorr_in_highCV_normal,'.--','LineWidth',2)
% plot(dorders_wtdDist,perc_sigCorr_in_highCV_wtdDist,'.-','LineWidth',2)
% plot(dorders_expDecay100,perc_sigCorr_in_highCV_expDecay100,'.--','LineWidth',2)
% plot(dorders_expDecay2,perc_sigCorr_in_highCV_expDecay2,'.--','LineWidth',2)
% plot(dorders_wtdDist_expDecay100,perc_sigCorr_in_highCV_wtdDist_expDecay100,'.-','LineWidth',2)
% plot(dorders_wtdDist_expDecay2,perc_sigCorr_in_highCV_wtdDist_expDecay2,'.-','LineWidth',2)
% xlabel('most global --> most local (scale not comparable!)');
% ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');
% legend({'oriDist + oriDecay','wtdDist + oriDecay',...
%     'oriDist + wtdDecay (100)','oriDist + wtdDecay (2)'...
%     'wtdDist + wtdDecay (100)','wtdDist + wtdDecay (2)'},'Location','north','FontSize',10,'FontWeight','bold');
%% plot the scatter plot of r of correlation vs. final improved FPA 
normalizedLevel = normalizedLevel_pro_perPro;
valid_rxns = valid_rxns_pro_perPro ;
%Computing the correlation between a reaction expression and measured growth rate
r_exp_only=[];
p_exp_only = [];
testedRxn_exp_only = {};
rxnLabel = fluxTbl.Model_Reaction_ID;
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn_exp_only(end+1) = rxnLabel(j);
        % [rho(end+1),p_rho(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
        [r_exp_only(end+1),p_exp_only(end+1)] = corr(normalizedLevel(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        % 
        %Correcting for multiple hypothesis using FDR and significance
        %level of?
    end
end
fdr_r_exp_only = mafdr(p_exp_only,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r_exp_only(fdr_r_exp_only<0.05)>0));


% correlation of FPA prediction
load(['output/Titration_relativeExp_wtdDist_expDecay.mat'])
nn = 9;
zz = 2;
FP = FP_collection_2{zz}{nn};
relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
for i = 1:size(FP,1)
    for j = 1:(size(FP,2)-1)
        if  mean(fluxMat(i,:)) > 0 
            if ~isnan(FP{i,j}(1))
                relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            elseif ~isnan(FP{i,j}(2))
                relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
            else
                fprintf('base = %f, n = %d, i = %d, j = %d, is a failed FPA\n',...
                    base(zz), dorders(nn),i,j);
                warning('check');
            end
        else
            if ~isnan(FP{i,j}(2))
                relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
            elseif ~isnan(FP{i,j}(1))
                relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            else
                error('check');
            end
        end
    end
end
rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
relFP = relFP(~rmInd,:);
valid_rxns = targetRxns(~rmInd);
%Computing the correlation between a reaction expression and measured growth rate
r_FPA=[];
deltaminmax = [];
testedRxn_FPA = {};
p_FPA = [];
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn_FPA(end+1) = rxnLabel(j);
        [r_FPA(end+1),p_FPA(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end
            
fdr_r_FPA = mafdr(p_FPA,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r_FPA(fdr_r_FPA<0.05)>0));

commonRxns = intersect(testedRxn_exp_only, testedRxn_FPA);
[A B] = ismember(commonRxns, testedRxn_exp_only);
r_exp_only = r_exp_only(B(A))';
p_exp_only = fdr_r_exp_only(B(A))';
[A B] = ismember(commonRxns, testedRxn_FPA);
r_FPA_extra = r_FPA(setdiff(1:length(r_FPA),B(A)))';
r_FPA_ori = r_FPA;
p_FPA_ori = fdr_r_FPA;
r_FPA = r_FPA(B(A))';
p_FPA = fdr_r_FPA(B(A))';

plot(r_exp_only,r_FPA,'.k','MarkerSize',10)
cutoff = 0.25;
hold on 
plot([-1 1], [-1 1],'--')
% plot(r_exp(r_exp>0.3 & r_exp - r_FPA > 0.3 ),r_FPA(r_exp>0.3 & (r_exp - r_FPA) > 0.3),'b.')
% plot(r_exp(r_FPA>0.3 & r_FPA - r_exp > 0.3 ),r_FPA(r_FPA>0.3 & (r_FPA - r_exp) > 0.3),'r.')
sig1 = p_exp_only<0.05 & r_exp_only >0;
sig2 = p_FPA<0.05 & r_FPA >0;
plot(r_exp_only(sig1 & ~sig2),r_FPA(sig1 & ~sig2),'b.','MarkerSize',10)
plot(r_exp_only(~sig1 & sig2),r_FPA(~sig1 & sig2),'r.','MarkerSize',10)
hold off
xlabel('PCC by relative expression');
ylabel('PCC by relative flux potential');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3, 2.7];
plt.LineWidth = 1;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/r_comparison_FPA_expressionOnly.tiff']);

histogram(r_FPA_extra,'FaceColor','none')
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.75, 1.5];
plt.LineWidth = 1;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.export(['figures/r_comparison_FPA_expressionOnly_extra_rxns_FPA.tiff']);

% x = ones(length(r_FPA_extra),1);
% y = r_FPA_extra;
% beeswarm(x,y,'sort_style','down')


%% target out -- local information's effect 
% correlation of FPA prediction
load(['output/Titration_relativeExp_wtdDist_expDecay_base2_n6_7_targetOut.mat'])
nn = 1;
FP = FP_collection_2{nn};
relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
for i = 1:size(FP,1)
    for j = 1:(size(FP,2)-1)
        if  mean(fluxMat(i,:)) > 0 
            if ~isnan(FP{i,j}(1))
                relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            elseif ~isnan(FP{i,j}(2))
                relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
            else
                fprintf('base = %f, n = %d, i = %d, j = %d, is a failed FPA\n',...
                    base(zz), dorders(nn),i,j);
                warning('check');
            end
        else
            if ~isnan(FP{i,j}(2))
                relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
            elseif ~isnan(FP{i,j}(1))
                relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
            else
                error('check');
            end
        end
    end
end
rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
relFP = relFP(~rmInd,:);
valid_rxns = targetRxns(~rmInd);
%Computing the correlation between a reaction expression and measured growth rate
r_FPA_targetOut=[];
deltaminmax = [];
testedRxn_FPA_targetOut = {};
p_FPA_targetOut = [];
for j = 1:length(rxnLabel)
    fluxMeasure = fluxMat_normalized(j,:);
    if any(strcmp(valid_rxns,rxnLabel{j}))
        testedRxn_FPA_targetOut(end+1) = rxnLabel(j);
        [r_FPA_targetOut(end+1),p_FPA_targetOut(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
        %Correcting for multiple hypothesis using FDR and significance level of
        deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
    end
end
            
fdr_r_FPA_targetOut = mafdr(p_FPA_targetOut,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(r_FPA_targetOut(fdr_r_FPA_targetOut<0.05)>0));

[A B] = ismember(commonRxns, testedRxn_FPA_targetOut);
r_FPA_targetOut_extra = r_FPA_targetOut(setdiff(1:length(r_FPA_targetOut),B(A)))';
r_FPA_targetOut_ori = r_FPA_targetOut;
p_FPA_targetOut_ori = fdr_r_FPA_targetOut;
r_FPA_targetOut = r_FPA_targetOut(B(A))';
p_FPA_targetOut = fdr_r_FPA_targetOut(B(A))';



plot(r_exp_only,r_FPA_targetOut,'.k','MarkerSize',10)
cutoff = 0.25;
hold on 
plot([-1 1], [-1 1],'--')
% plot(r_exp(r_exp>0.3 & r_exp - r_FPA > 0.3 ),r_FPA(r_exp>0.3 & (r_exp - r_FPA) > 0.3),'b.')
% plot(r_exp(r_FPA>0.3 & r_FPA - r_exp > 0.3 ),r_FPA(r_FPA>0.3 & (r_FPA - r_exp) > 0.3),'r.')
sig1 = p_exp_only<0.05 & r_exp_only >0;
sig2 = p_FPA_targetOut<0.05 & r_FPA_targetOut >0;
plot(r_exp_only(sig1 & ~sig2),r_FPA_targetOut(sig1 & ~sig2),'b.','MarkerSize',10)
plot(r_exp_only(sig2 & ~sig1),r_FPA_targetOut(sig2 & ~sig1),'r.','MarkerSize',10)
hold off
xlabel('PCC by relative expression');
ylabel('PCC by relative flux potential');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3, 2.7];
plt.LineWidth = 1;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/r_comparison_FPAtargetOut_expressionOnly.tiff']);

histogram(r_FPA_targetOut_extra,'FaceColor','none')
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.75, 1.5];
plt.LineWidth = 1;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.export(['figures/r_comparison_FPAtargetOut_expressionOnly_extra_rxns_FPA.tiff']);


plot(r_FPA_ori,r_FPA_targetOut_ori,'.k','MarkerSize',10)
cutoff = 0.25;
hold on 
plot([-1 1], [-1 1],'--')
% plot(r_exp(r_exp>0.3 & r_exp - r_FPA > 0.3 ),r_FPA(r_exp>0.3 & (r_exp - r_FPA) > 0.3),'b.')
% plot(r_exp(r_FPA>0.3 & r_FPA - r_exp > 0.3 ),r_FPA(r_FPA>0.3 & (r_FPA - r_exp) > 0.3),'r.')
sig1 = p_FPA_ori<0.05 & r_FPA_ori >0;
sig2 = p_FPA_targetOut_ori<0.05 & r_FPA_targetOut_ori >0;
plot(r_FPA_ori(sig1 & ~sig2),r_FPA_targetOut_ori(sig1 & ~sig2),'b.','MarkerSize',10)
plot(r_FPA_ori(sig2 & ~sig1),r_FPA_targetOut_ori(sig2 & ~sig1),'r.','MarkerSize',10)
hold off
xlabel('PCC by target in');
ylabel('PCC by target out');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3, 2.7];
plt.LineWidth = 1;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.export(['figures/r_comparison_FPAtargetOut_FPA.tiff']);
%% ######################################
%% ######################################
%% ######################################
% %% same comparsion by raw-raw data
% %% evaluate by the data cleanning method
% load Titration_rawExp_oriDist_oriDecay.mat
% %%
% load Titration_rawExp_wtdDist_oriDecay.mat
% %% correlation for 1-d titration
% dorders = n;
% rhoMat = zeros(length(targetRxns),length(dorders));
% N_sigCorr = zeros(length(dorders),1);
% N_sigCorr_highCV = zeros(length(dorders),1);
% perc_sigCorr_in_highCV = zeros(length(dorders),1);
% for nn = 1: length(dorders)
% %%
% % nn = 10;
% dorders(nn)
% FP = FP_collection{nn};
% relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
% for i = 1:size(FP,1)
%     for j = 1:(size(FP,2)-1)
%         if  mean(fluxMat(i,:)) > 0 
%             if ~isnan(FP{i,j}(1))
%                 relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
%             elseif ~isnan(FP{i,j}(2))
%                 relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
%             else
%                 error('check');
%             end
%         else
%             if ~isnan(FP{i,j}(2))
%                 relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
%             elseif ~isnan(FP{i,j}(1))
%                 relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
%             else
%                 error('check');
%             end
%         end
%     end
% end
% rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
% relFP = relFP(~rmInd,:);
% valid_rxns = targetRxns(~rmInd);
% %Computing the correlation between a reaction expression and measured growth rate
% rho=[];
% p=[];
% cv = [];
% testedRxn = {};
% for j = 1:length(rxnLabel)
%     fluxMeasure = fluxMat(j,:);
%     if any(strcmp(valid_rxns,rxnLabel{j}))
%         testedRxn(end+1) = rxnLabel(j);
%         [rho(end+1),p(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
%         %Correcting for multiple hypothesis using FDR and significance level of
%         cv(end+1) = std(relFP(strcmp(valid_rxns,rxnLabel{j}),:))/mean(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
%     end
% end
% sum(rho(p<0.05)>0)
% sum(rho(p<0.05)>0)/length(testedRxn)
% histogram(rho,'BinLimits',[-1,1],'NumBins',20)
% [A B] = ismember(testedRxn,targetRxns);
% rhoMat(B(A),nn) = rho;
% N_sigCorr(nn) = sum(rho(p<0.05)>0);
% N_sigCorr_highCV(nn) = sum(rho(p<0.05 & cv > 0.1)>0);
% perc_sigCorr_in_highCV(nn) = sum(rho(p<0.05 & cv > 0.1)>0) / sum(cv > 0.1);
% end
% %% analyze some overall metric
% figure(1)
% plot(dorders,N_sigCorr,'.-')
% xlabel('distance order');
% ylabel('number of sig corr rxns');
% figure(2)
% plot(dorders,N_sigCorr_highCV,'.-')
% xlabel('distance order');
% ylabel('number of sig-corr & rFP-CV > 0.1 rxns ');
% figure(3)
% plot(dorders,perc_sigCorr_in_highCV,'.-')
% xlabel('distance order');
% ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');
% 
% 
% %%
% load Titration_rawExp_oriDist_expDecay.mat
% %%
% load Titration_rawExp_wtdDist_expDecay.mat
% %% correlation for 2-d titration
% N_sigCorr = zeros(length(n2),length(base));
% N_sigCorr_highCV = zeros(length(n2),length(base));
% perc_sigCorr_in_highCV = zeros(length(n2),length(base));
% for zz = 1:length(base)
%     dorders = n2;
%     rhoMat = zeros(length(targetRxns),length(dorders));
%     for nn = 1: length(dorders)
%     %%
%     % nn = 10;
%     %dorders(nn)
%     FP = FP_collection_2{zz}{nn};
%     relFP = nan(size(FP,1),size(FP,2)-1);%we choose one from f and r as a prediction
%     for i = 1:size(FP,1)
%         for j = 1:(size(FP,2)-1)
%             if  mean(fluxMat(i,:)) > 0 
%                 if ~isnan(FP{i,j}(1))
%                     relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
%                 elseif ~isnan(FP{i,j}(2))
%                     relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
%                 else
%                     fprintf('base = %f, n = %d, i = %d, j = %d, is a failed FPA\n',...
%                         base(zz), dorders(nn),i,j);
%                     warning('check');
%                 end
%             else
%                 if ~isnan(FP{i,j}(2))
%                     relFP(i,j) = FP{i,j}(2) ./ FP{i,end}(2);
%                 elseif ~isnan(FP{i,j}(1))
%                     relFP(i,j) = FP{i,j}(1) ./ FP{i,end}(1);
%                 else
%                     error('check');
%                 end
%             end
%         end
%     end
%     rmInd = all(isnan(relFP),2) | all(relFP == relFP(:,1),2);
%     relFP = relFP(~rmInd,:);
%     valid_rxns = targetRxns(~rmInd);
%     %Computing the correlation between a reaction expression and measured growth rate
%     rho=[];
%     p=[];
%     cv = [];
%     testedRxn = {};
%     for j = 1:length(rxnLabel)
%         fluxMeasure = fluxMat(j,:);
%         if any(strcmp(valid_rxns,rxnLabel{j}))
%             testedRxn(end+1) = rxnLabel(j);
%             [rho(end+1),p(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Spearman');
%             %Correcting for multiple hypothesis using FDR and significance level of
%             cv(end+1) = std(relFP(strcmp(valid_rxns,rxnLabel{j}),:))/mean(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
%         end
%     end
%     sum(rho(p<0.05)>0);
%     sum(rho(p<0.05)>0)/length(testedRxn);
%     histogram(rho,'BinLimits',[-1,1],'NumBins',20);
%     [A B] = ismember(testedRxn,targetRxns);
%     rhoMat(B(A),nn) = rho;
%     N_sigCorr(nn,zz) = sum(rho(p<0.05)>0);
%     N_sigCorr_highCV(nn,zz) = sum(rho(p<0.05 & cv > 0.1)>0);
%     perc_sigCorr_in_highCV(nn,zz) = sum(rho(p<0.05 & cv > 0.1)>0) / sum(cv > 0.1);
%     end
% end
% %% analyze some overall metric
% figure(4)
% hold on 
% for zz = 1:length(base)
%     plot(dorders,N_sigCorr(:,zz),'.-')
% end
% hold off
% xlabel('distance order');
% ylabel('number of sig corr rxns');
% legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))
% 
% figure(5)
% hold on 
% for zz = 1:length(base)
%     plot(dorders,N_sigCorr_highCV(:,zz),'.-')
% end
% hold off
% xlabel('distance order');
% ylabel('number of sig-corr & rFP-CV > 0.1 rxns ');
% legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))
% 
% figure(6)
% hold on 
% for zz = 1:length(base)
%     plot(dorders,perc_sigCorr_in_highCV(:,zz),'.-')
% end
% hold off
% xlabel('distance order');
% ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');
% legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))
% 
% 
% figure(7)
% hold on 
% for zz = 1:length(base)
%     plot(dorders(1:13),N_sigCorr_highCV(1:13,zz),'.-')
% end
% hold off
% xlabel('distance order');
% ylabel('number of sig-corr & rFP-CV > 0.1 rxns ');
% legend(cellfun(@(x) ['base = ',x],strsplit(num2str(base)), 'UniformOutput',0))
% 
% %% plot four setups together: 1-d tritrations
% dorders_normal = dorders/max(dorders);
% N_sigCorr_normal = N_sigCorr;
% N_sigCorr_highCV_normal = N_sigCorr_highCV;
% perc_sigCorr_in_highCV_normal = perc_sigCorr_in_highCV;
% 
% dorders_wtdDist = dorders/max(dorders);
% N_sigCorr_wtdDist = N_sigCorr;
% N_sigCorr_highCV_wtdDist = N_sigCorr_highCV;
% perc_sigCorr_in_highCV_wtdDist = perc_sigCorr_in_highCV;
% %% exp decay only -2d
% dorders_expDecay100 = 1-dorders/max(dorders);
% N_sigCorr_expDecay100 = N_sigCorr(:,7);% base =100
% N_sigCorr_highCV_expDecay100 = N_sigCorr_highCV(:,7);
% perc_sigCorr_in_highCV_expDecay100 = perc_sigCorr_in_highCV(:,7);
% 
% dorders_expDecay2 = 1-dorders/max(dorders);
% N_sigCorr_expDecay2 = N_sigCorr(:,2);% base =2
% N_sigCorr_highCV_expDecay2 = N_sigCorr_highCV(:,2);
% perc_sigCorr_in_highCV_expDecay2 = perc_sigCorr_in_highCV(:,2);
% %% exp decay + wtd dist -2d
% dorders_wtdDist_expDecay100 = 1-dorders/max(dorders);
% N_sigCorr_wtdDist_expDecay100 = N_sigCorr(:,7);% base =100
% N_sigCorr_highCV_wtdDist_expDecay100 = N_sigCorr_highCV(:,7);
% perc_sigCorr_in_highCV_wtdDist_expDecay100 = perc_sigCorr_in_highCV(:,7);
% 
% dorders_wtdDist_expDecay2 = 1-dorders/max(dorders);
% N_sigCorr_wtdDist_expDecay2 = N_sigCorr(:,2);% base =2
% N_sigCorr_highCV_wtdDist_expDecay2 = N_sigCorr_highCV(:,2);
% perc_sigCorr_in_highCV_wtdDist_expDecay2 = perc_sigCorr_in_highCV(:,2);
% %% plot 
% figure(4)
% hold on
% plot(dorders_normal,N_sigCorr_normal,'.--','LineWidth',2)
% plot(dorders_wtdDist,N_sigCorr_wtdDist,'.-','LineWidth',2)
% plot(dorders_expDecay100,N_sigCorr_expDecay100,'.--','LineWidth',2)
% plot(dorders_expDecay2,N_sigCorr_expDecay2,'.--','LineWidth',2)
% plot(dorders_wtdDist_expDecay100,N_sigCorr_wtdDist_expDecay100,'.-','LineWidth',2)
% plot(dorders_wtdDist_expDecay2,N_sigCorr_wtdDist_expDecay2,'.-','LineWidth',2)
% xlabel('most global --> most local (scale not comparable!)');
% ylabel('number of sig corr rxns');
% legend({'oriDist + oriDecay','wtdDist + oriDecay',...
%     'oriDist + wtdDecay (100)','oriDist + wtdDecay (2)'...
%     'wtdDist + wtdDecay (100)','wtdDist + wtdDecay (2)'},'Location','north','FontSize',10,'FontWeight','bold');
% 
% figure(5)
% hold on
% plot(dorders_normal,N_sigCorr_highCV_normal,'.--','LineWidth',2)
% plot(dorders_wtdDist,N_sigCorr_highCV_wtdDist,'.-','LineWidth',2)
% plot(dorders_expDecay100,N_sigCorr_highCV_expDecay100,'.--','LineWidth',2)
% plot(dorders_expDecay2,N_sigCorr_highCV_expDecay2,'.--','LineWidth',2)
% plot(dorders_wtdDist_expDecay100,N_sigCorr_highCV_wtdDist_expDecay100,'.-','LineWidth',2)
% plot(dorders_wtdDist_expDecay2,N_sigCorr_highCV_wtdDist_expDecay2,'.-','LineWidth',2)
% xlabel('most global --> most local (scale not comparable!)');
% ylabel('number of sig-corr & rFP-CV > 0.1 rxns ');
% legend({'oriDist + oriDecay','wtdDist + oriDecay',...
%     'oriDist + wtdDecay (100)','oriDist + wtdDecay (2)'...
%     'wtdDist + wtdDecay (100)','wtdDist + wtdDecay (2)'},'Location','north','FontSize',10,'FontWeight','bold');
% 
% figure(6)
% hold on
% plot(dorders_normal,perc_sigCorr_in_highCV_normal,'.--','LineWidth',2)
% plot(dorders_wtdDist,perc_sigCorr_in_highCV_wtdDist,'.-','LineWidth',2)
% plot(dorders_expDecay100,perc_sigCorr_in_highCV_expDecay100,'.--','LineWidth',2)
% plot(dorders_expDecay2,perc_sigCorr_in_highCV_expDecay2,'.--','LineWidth',2)
% plot(dorders_wtdDist_expDecay100,perc_sigCorr_in_highCV_wtdDist_expDecay100,'.-','LineWidth',2)
% plot(dorders_wtdDist_expDecay2,perc_sigCorr_in_highCV_wtdDist_expDecay2,'.-','LineWidth',2)
% xlabel('most global --> most local (scale not comparable!)');
% ylabel('percentage of sig-corr rxns in all reaction with rFP-CV > 0.1');
% legend({'oriDist + oriDecay','wtdDist + oriDecay',...
%     'oriDist + wtdDecay (100)','oriDist + wtdDecay (2)'...
%     'wtdDist + wtdDecay (100)','wtdDist + wtdDecay (2)'},'Location','north','FontSize',10,'FontWeight','bold');
