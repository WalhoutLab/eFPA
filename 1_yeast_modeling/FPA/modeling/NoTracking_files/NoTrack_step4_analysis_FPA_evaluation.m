%% About
% FPA evaluation. systematically compare the FPA prediction with no
% integration and evalute the distance boundaries
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/

%% 1. load the model and prepare the model
addpath('./../scripts/')
model = loadYeatModel();
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
%% 3. load the expression files, distance matrix, and other optional inputs
% load the distance matrix
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
extRxns = model.rxns(findExcRxns(model));
manualPenalty = extRxns;
manualPenalty(:,2) = mat2cell(zeros(length(manualPenalty),1),ones(length(manualPenalty),1));
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
fluxMat_normalized = fluxMat;
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
%% load expression levels
load('output/normalizedLevels_partialExcluded.mat');
%% precalc penalty 
penalty_pro = ones(length(model.rxns),size(fluxMat,2)+1);
[A B] = ismember(model.rxns,valid_rxns_pro_perPro);
penalty_pro(A,1:(end-1)) = ones(size(normalizedLevel_pro_perPro(B(A),:),1),size(normalizedLevel_pro_perPro(B(A),:),2)) ./ normalizedLevel_pro_perPro(B(A),:);
% apply additional penalty to the exchange of limiting nutrients
for i = 1:size(manualPenalty,1)
    penalty_pro(strcmp(model.rxns,manualPenalty{i,1}),:) = manualPenalty{i,2};
end
penalty_pro(strcmp(model.rxns,'r_2005'),1:5) = 10;
penalty_pro(strcmp(model.rxns,'r_1714'),6:10) = 10;
penalty_pro(strcmp(model.rxns,'r_1654'),11:15) = 10;
%%
minmaxCutoff = 0.1;
mySet = intersect(valid_rxns_pro_perPro, rxnLabel);
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
        dorders(nn);
        FP = FP_collection{nn};
        relFP = nan(size(FP,1),size(FP,2)-1);% we choose one from f and r as a prediction
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
        % Computing the correlation between a reaction expression and measured growth rate
        r=[];
        p_r=[];
        deltaminmax = [];
        testedRxn = {};
        for j = 1:length(rxnLabel)
            fluxMeasure = fluxMat_normalized(j,:);
            if any(strcmp(valid_rxns,rxnLabel{j}))
                testedRxn(end+1) = rxnLabel(j);
                [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
                % to avoid false positives related to numeric issues, we
                % filter by the range of the rFP values
                deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
            end
        end

        fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
        fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & deltaminmax > minmaxCutoff & ismember(testedRxn, mySet))>0));

        if (dorders(nn) == 1.5)
            figure;
            hold on
            histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
            xlim([-1,1]);
            xlabel('Correlation coefficient');
            ylabel('Number of reactions');
            histogram(r(fdr_r<0.05 & r >0 & deltaminmax > minmaxCutoff),'FaceColor','#D95319','BinEdges',-1:0.2:1)
            legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
            hold off
            plt = Plot(); % create a Plot object and grab the current figure
            plt.BoxDim = [1.95, 1.6125];
            plt.LineWidth = 1;
            plt.FontSize = 10;
            plt.XTick = -1:0.2:1;
            plt.XMinorTick = 'off';
            plt.YMinorTick = 'off';
            plt.ShowBox = 'off';
            plt.LegendLoc = 'NorthWest';
            plt.FontName = 'Arial';
        end

        [A B] = ismember(testedRxn,targetRxns);
        rMat(B(A),nn) = r;
        N_sigCorr(nn) = sum(r(fdr_r<0.05)>0);
        N_sigCorr_highCV(nn) = sum(fdr_r<0.05 & deltaminmax > minmaxCutoff & r > 0 & ismember(testedRxn, mySet));
        perc_sigCorr_in_highCV(nn) = sum(r(fdr_r<0.05 & deltaminmax > minmaxCutoff)>0) / sum(deltaminmax > minmaxCutoff);
    end
    %% save data to plot four setups together: 1-d tritrations
    if zz==1
        dorders_normal = dorders;% dorders/max(dorders);
        N_sigCorr_normal = N_sigCorr;
        N_sigCorr_highCV_normal = N_sigCorr_highCV;
        perc_sigCorr_in_highCV_normal = perc_sigCorr_in_highCV;
    elseif zz ==2
        dorders_wtdDist = dorders;% dorders/max(dorders);
        N_sigCorr_wtdDist = N_sigCorr;
        N_sigCorr_highCV_wtdDist = N_sigCorr_highCV;
        perc_sigCorr_in_highCV_wtdDist = perc_sigCorr_in_highCV;
    end
end

% the 2-d titrition data
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
            % Computing the correlation
            r=[];
            p_r=[];
            deltaminmax = [];
            testedRxn = {};
            for j = 1:length(rxnLabel)
                fluxMeasure = fluxMat_normalized(j,:);
                if any(strcmp(valid_rxns,rxnLabel{j}))
                    testedRxn(end+1) = rxnLabel(j);
                    [r(end+1),p_r(end+1)] = corr(relFP(strcmp(valid_rxns,rxnLabel{j}),:)',abs(fluxMeasure)','type','Pearson');
                    deltaminmax(end+1) = max(relFP(strcmp(valid_rxns,rxnLabel{j}),:)) - min(relFP(strcmp(valid_rxns,rxnLabel{j}),:));
                end
            end
            
            fdr_r = mafdr(p_r,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
            fprintf('%d rxns give significant positive correlation by pearson\n',sum(r(fdr_r<0.05 & deltaminmax > minmaxCutoff & ismember(testedRxn, mySet))>0));
      
            if (dorders(nn) == 6 && base(zz) ==2)
                figure;
                hold on
                histogram(r,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
                xlim([-1,1]);
                xlabel('Correlation coefficient');
                ylabel('Number of reactions');
                histogram(r(fdr_r<0.05 & r >0 & deltaminmax > minmaxCutoff),'FaceColor','#D95319','BinEdges',-1:0.2:1)
                legend({'all testable reactions',sprintf('significantly correlated \nreactions')})
                hold off
                plt = Plot(); % create a Plot object and grab the current figure
                plt.BoxDim = [1.95, 1.6125];
                plt.LineWidth = 1;
                plt.FontSize = 10;
                plt.XTick = -1:0.2:1;
                plt.XMinorTick = 'off';
                plt.YMinorTick = 'off';
                plt.ShowBox = 'off';
                plt.LegendLoc = 'NorthWest';
                plt.FontName = 'Arial';
            end

            [A B] = ismember(testedRxn,targetRxns);
            rMat(B(A),nn) = r;
            N_sigCorr(nn,zz) = sum(r(fdr_r<0.05)>0);
            N_sigCorr_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > minmaxCutoff & ismember(testedRxn, mySet))>0);
            perc_sigCorr_in_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > minmaxCutoff)>0) / sum(deltaminmax > minmaxCutoff);
        end
    end
    %% analyze some overall metric
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
    plt.BoxDim = [3, 2.4];
    plt.LineWidth = 1;
    plt.FontSize = 10;
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.ShowBox = 'off';
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    
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
    plt.BoxDim = [3, 2.4];
    plt.LineWidth = 1;
    plt.FontSize = 10;
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.ShowBox = 'off';
    plt.FontName = 'Arial';
    if jj ==1
        plt.LegendLoc = 'NorthEast';
    elseif jj ==2
        plt.LegendLoc = 'South';
    end
    plt.FontName = 'Arial';
    
     if jj==1
        %% exp decay only -2d
        dorders_expDecay1000 = dorders;
        N_sigCorr_expDecay1000 = N_sigCorr(:,8);% base =1000
        N_sigCorr_highCV_expDecay1000 = N_sigCorr_highCV(:,8);
        perc_sigCorr_in_highCV_expDecay1000 = perc_sigCorr_in_highCV(:,8);

        dorders_expDecay2 = dorders;
        N_sigCorr_expDecay2 = N_sigCorr(:,2);% base =2
        N_sigCorr_highCV_expDecay2 = N_sigCorr_highCV(:,2);
        perc_sigCorr_in_highCV_expDecay2 = perc_sigCorr_in_highCV(:,2);
     elseif jj ==2
    %% exp decay + wtd dist -2d
        dorders_wtdDist_expDecay1000 = dorders;% 
        N_sigCorr_wtdDist_expDecay1000 = N_sigCorr(:,8);% base =1000
        N_sigCorr_highCV_wtdDist_expDecay1000 = N_sigCorr_highCV(:,8);
        perc_sigCorr_in_highCV_wtdDist_expDecay1000 = perc_sigCorr_in_highCV(:,8);

        dorders_wtdDist_expDecay2 = dorders;% 
        N_sigCorr_wtdDist_expDecay2 = N_sigCorr(:,2);% base =2
        N_sigCorr_highCV_wtdDist_expDecay2 = N_sigCorr_highCV(:,2);
        perc_sigCorr_in_highCV_wtdDist_expDecay2 = perc_sigCorr_in_highCV(:,2);
     end
end

%% plot 
baseline = 46; % number of correlated rxns by ROI expression only

figure;
hold on
plot(dorders_wtdDist_expDecay2(1:14),N_sigCorr_highCV_wtdDist_expDecay2(1:14) ,'o-','LineWidth',2,'Color','#D95319','MarkerSize', 3)
plot(dorders_wtdDist_expDecay1000(1:14),N_sigCorr_highCV_wtdDist_expDecay1000(1:14) ,'o-','LineWidth',2,'Color','#EDB120','MarkerSize', 3)
plot(dorders_wtdDist(1:18),N_sigCorr_highCV_wtdDist(1:18) ,'o-','LineWidth',2,'Color','#0072BD','MarkerSize', 3)
plot(dorders_normal(1:18),N_sigCorr_highCV_normal(1:18) ,'o-k','LineWidth',2,'MarkerSize', 3)
yline(baseline,'--','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Distance order/boundary');
ylabel('Number of significantly correlated reactions ');
ylim([20 55])
legend({'improved FPA (base = 2)','improved FPA (base = 1000)','original FPA + weighted distance', 'original FPA','expression only'},'FontSize',7);
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2.85, 2.35];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.LegendLoc = 'SouthEast';
plt.export(['figures/NoTrack_FPAbenchmarking_156_',num2str(minmaxCutoff),'.pdf']);