%%
initCobraToolbox
%%
addpath ./../input/YeastJoshua/
addpath ./../scripts/
addpath ./../../bins/
addpath('./../scripts/')
addpath('./../../../PlotPub/lib/')
addpath(genpath('./../other_methods/REMI/'));

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


proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);

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

%% evaluate different FPA
% load result of original FPA
setups = {'oriDist_oriDecay'};
for zz = 1
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
            plt.BoxDim = [1.95, 1.6125];
            plt.LineWidth = 1;
            plt.FontSize = 10;
            plt.XTick = -1:0.2:1;
            plt.XMinorTick = 'off';
            plt.YMinorTick = 'off';
            plt.ShowBox = 'off';
            plt.LegendLoc = 'NorthWest';
            plt.FontName = 'Arial';

            n_corr_ori_FPA = sum(fdr_r<0.05 & r >0 & deltaminmax > 0.2);
        end

         if (dorders(nn) == 0) % negative control showing integrating everything
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
            plt.BoxDim = [1.95, 1.6125];
            plt.LineWidth = 1;
            plt.FontSize = 10;
            plt.XTick = -1:0.2:1;
            plt.XMinorTick = 'off';
            plt.YMinorTick = 'off';
            plt.ShowBox = 'off';
            plt.LegendLoc = 'NorthWest';
            plt.FontName = 'Arial';

            n_corr_no_distance_decay = sum(fdr_r<0.05 & r >0 & deltaminmax > 0.2);
        end

        [A B] = ismember(testedRxn,targetRxns);
        rMat(B(A),nn) = r;
        N_sigCorr(nn) = sum(r(fdr_r<0.05)>0);
        N_sigCorr_highCV(nn) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0);
        perc_sigCorr_in_highCV(nn) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0) / sum(deltaminmax > 0.2);
    end
end

% load result of eFPA
setups = {'wtdDist_expDecay'};
for jj = 1

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
            plt.BoxDim = [1.95, 1.6125];
            plt.LineWidth = 1;
            plt.FontSize = 10;
            plt.XTick = -1:0.2:1;
            plt.XMinorTick = 'off';
            plt.YMinorTick = 'off';
            plt.ShowBox = 'off';
            plt.LegendLoc = 'NorthWest';
            plt.FontName = 'Arial';
        
            n_corr_eFPA_local = sum(fdr_r<0.05 & r >0 & deltaminmax > 0.2);
        end
        
        [A B] = ismember(testedRxn,targetRxns);
        rMat(B(A),nn) = r;
        N_sigCorr(nn,zz) = sum(r(fdr_r<0.05)>0);
        N_sigCorr_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0);
        perc_sigCorr_in_highCV(nn,zz) = sum(r(fdr_r<0.05 & deltaminmax > 0.2)>0) / sum(deltaminmax > 0.2);
    end
end

end

%% make the REMI matrix 
% will be seperated into a function

% We run REMI with its own GRP parsing algorithm and just supply the
% fold-changes

% ExpTable = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)

% we optimize the REMI parameter for fold-change thresholding in case it is
% important 
% cutoff    ncorr
% cutoff    ncorr
% 1         17
% 1.1       28
% 1.2       47
% 1.3       32
% 1.5       25
% 1.6       2
% 1.7       3
% 1.8       8
% 1.9       6
% 2         2

cutoffs = 1.2;
n_corr_REMI = test_REMI(cutoffs);


%% compass integration 

% we run compass in a custom manner because it cannot be used directly
% to be comparable as much as possible, we directly used the reaction
% expression level from FPA; the penalty is the inverse of it (1/x instead of
% 1/(x+1)). This should still be conceptually the same as compass
% calculation. other main algorithm was implemented in house 
 
% load compass result 
% load('/Users/xuhangli/Library/CloudStorage/OneDrive-UMassChanMedicalSchool/Walhout_Lab/eFPA_project/improvedFPA_final_run/1_yeast_modeling/FPA/modeling/COMPASS_result/Titration_relativeExp_oriDist_COMPASSLIKE.mat')

% the first one (n=0) is a compass setup simulation 
% we adopted the log transformation of the resistance score in compass
% (super condition is not used - we dont do the ratio here because minus
% log makes more sense)
% we ignore the meta reaction processing

parpool(8);

testN = [0 2.5]; % use 2.5 as a relative local pathway when using original decay 

n_corr_list = test_Compass(testN);
n_corr_COMPASS = n_corr_list(1);
n_corr_COMPASS_local_decay = n_corr_list(2);

%% delta FBA

% we did an honest try on the parameter optimizations, including cut_off
% (0.05, 0.1, 0.25) and constraints (w/ biomass, w/o biomass, w/ default
% exchange constriants, with free exchange constraints), maxflux_val (
% 10, 200) and epsilon (0.1, 10) and cannot get any good predictions. the
% correlation is always very poor. we ended up with the parameters similar
% to the showcase of muscle model in the original publication

% in theroy, we should use statement B for correlation analysis. however,
% it is quite intensive so we can only try on a few parameter choices. The
% best case of statement A is ~10 reactions predicted (with biomass, unlimited exchange
% and statement A, 5% cutoff). But this one is very difficult to reproduce
% because the many solving was cutoff by runtime and correlated reactions
% looks more like by chance. 

% so we ended up with a random setup since they are basically similar; we
% aligned the model constraints with eFPA (mostly unconstrained) and 
% we used 0.25 cutoff, vmax 200, espsilon 10 and statement A with MILP obj, 
% this setting is generally aligned with the setting used for muscle in the 
% original publication. 


% load('output/other_methods_deltaFBA_5percent.mat')
% load('output/other_methods_deltaFBA_25percent.mat')

n_corr_deltaFBA =  test_deltaFBA();

%%
baseline = 46;

% n_corr_REMI = 47;
% n_corr_deltaFBA = 12;

figure;
bar([n_corr_deltaFBA, n_corr_COMPASS, n_corr_no_distance_decay, n_corr_REMI, n_corr_COMPASS_local_decay,  n_corr_ori_FPA, n_corr_eFPA_local]);
title('method benchmarking');
xlabel('Method');
ylabel('Number of sig. corr. reactions');
set(gca, 'xticklabel', {'delta FBA', 'COMPASS', 'FPA - no distance decay', 'REMI',  'COMPASS + distance decay',  'original FPA', 'local-pathway eFPA'})
hold on; % This command allows you to add to the current plot
line(xlim, [46 46], 'Color', '#808080', 'LineStyle', '--')
hold off;
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


plt.export(['figures/FPAbenchmarking_agaist_others.pdf']);

%% to delete

% n_corr_REMI = [];
% model_ori = model;
% for zz = 1:length(cutoffs)
%     cutoff = cutoffs(zz)
%     model = model_ori;
% 
%     up_cutoff=cutoff; % more than 2 fold up regulated is known as up regulated genes (This paramter can change)
%     down_cutoff=1/cutoff;% less than 1/2 fold  is known as down regulated genes 
%     % now we are build gex Model:  which integeartes relative gene expression
%     refExp = ExpTable.P0_05; % random one 
%     REMImat = zeros(length(model.rxns),length(conditions));
%     %impose biomass
%     model = changeRxnBounds(model,'r_2111',0.1, 'l');
% 
%     for index = 1:length(conditions)
%         name = conditions{index};
%         myExp = ExpTable.(name);
%         fc = 2.^myExp;
%         [rxns,ratio]=evaluateGPR(model,ExpTable.Gene,fc,@geomean,@mean);
%         % find up- and down- regulated reactions
%         indUP=find(ratio>up_cutoff);
%         ratioUP=ratio(indUP);
%         indDOWN=find(ratio<down_cutoff);
%         ratioDOWN=ratio(indDOWN);
%         regRxnRatio=[ratioUP;ratioDOWN];
%         regRxns=[rxns(indUP)';rxns(indDOWN)'];
% 
%         % avoid numerical error (more than 100 fold is taken only 100)
%         regRxnRatio(regRxnRatio>100)=100;
%         regRxnRatio(regRxnRatio<1/100)=1/100;
% 
%         % if we want to add relative constraint into TFA model then we need to add
%         % net flux variable to the model using 'addNetFluxVariablesNEW'
%         % and if one want to add into FBA model  then evalute scripts given below
%         mTmp=model;
%         mTmp.rev=ones(numel(mTmp.rxns),1);
%         use_m=addUseVariablesDH(mTmp);
%         netM=addNetFluxVariablesNEW(use_m);
%         [gex]=addRelConsExpression(netM,netM,regRxns,regRxnRatio);
%         sol_gex=solveTFBAmodel(gex,false,'gurobi',true);
% 
%         remiOut = sol_gex;
%         [A B] = ismember(cellfun(@(x) ['PERTURB_NF_',x],model.rxns,'UniformOutput',false),gex.varNames);
%         [C D] = ismember(cellfun(@(x) ['NF_',x],model.rxns,'UniformOutput',false),gex.varNames);
%         REMImat(1:length(model.rxns),index) = abs(remiOut.x.x(B(A))) ./ abs(remiOut.x.x(D(C)));
% %         [A B] = ismember(cellfun(@(x) ['PERTURB_R_',x],model.rxns,'UniformOutput',false),gex.varNames);
% %         [C D] = ismember(cellfun(@(x) ['R_',x],model.rxns,'UniformOutput',false),gex.varNames);
% %         REMImat((1+length(model.rxns)):2*length(model.rxns),index) = remiOut.x.x(B(A)) ./ abs(remiOut.x.x(D(C)));
%     end
%     REMImat(isnan(REMImat)) = 1;
%     REMImat(isinf(REMImat)) = 100;
%     REMImat(REMImat==0) = 1e-2;
%     %% correlation
% 
%     fluxMat = fluxMat_normalized;
%     %Computing the correlation between a reaction expression and measured growth rate
%     rho=[];
%     p=[];
%     testedRxn = {};
%     rxnLabel = fluxTbl.Model_Reaction_ID;
%     valid_rxns = model.rxns;
%     for j = 1:length(rxnLabel)
%         fluxMeasure = fluxMat(j,:);
%         if any(strcmp(valid_rxns,rxnLabel{j}))
%             testedRxn(end+1) = rxnLabel(j);
%             prediction = REMImat(strcmp(valid_rxns,rxnLabel{j}),:)';
% 
% %             if mean(fluxMeasure) > 0
% %                 prediction = REMImat(strcmp(valid_rxns,rxnLabel{j}),:)';
% %             else
% %                 prediction = REMImat(find(strcmp(valid_rxns,rxnLabel{j})) + length(valid_rxns),:)';
% %             end
% 
%             [rho(end+1),p(end+1)] = corr(prediction,abs(fluxMeasure)','type','Pearson');
%             % Correcting for multiple hypothesis using FDR and significance level of
%         end
%     end
%     testedRxn_pro = testedRxn;
%     p_pro = p;
%     rho_pro = rho;
%     figure
%     histogram(rho)
%     sum(rho(p<0.05)>0)
% 
%     fdr = mafdr(p,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
%     n_corr_REMI(zz) = sum(fdr<0.05 & rho > 0);
% end
