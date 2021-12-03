%%
% PS. if we decouple ATP from the two reaction, the predictive power is
% gone 
%% load data
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
rxnLabel = fluxTbl.Model_Reaction_ID;
fluxMat_normalized = fluxMat;

% dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
% factor = repmat(dwTbl.gDCW_mL',size(fluxMat,1),1);
% fluxMat_normalized = fluxMat_normalized * 1000 ./ factor; %mmoles/hr/gDW

GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
model = loadYeatModel();
%% first extract the contributing rxns
load('output/Titration_relativeExp_wtdDist_expDecay_FineGrained.mat');
targetRxns = {'r_0460', 'r_0485'};
optimalBound = [18,18.5];
for kk = 1:length(targetRxns)
    targetRxn = targetRxns{kk};
    %% load the FPA table set 
    load(['fine_titration/base_1000_n_', num2str(optimalBound(kk)),'.mat']);
    %% clean up the predictions (only keep the relavant directions)
    FPsol = cell(size(FluxPotential_solutions,1),size(FluxPotential_solutions,2));%we choose one from f and r as a prediction
    for i = 1:size(FluxPotential_solutions,1)
        for j = 1:(size(FluxPotential_solutions,2))
            if  mean(fluxMat(i,:)) > 0 
                if ~isempty(FluxPotential_solutions{i,j}(1))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(1);
                elseif ~isempty(FluxPotential_solutions{i,j}(2))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(2);
                else
                    fprintf('i = %d, j = %d, is a failed FPA\n',...
                        i,j);
                    warning('check');
                end
            else
                if ~isempty(FluxPotential_solutions{i,j}(2))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(2);
                elseif ~isempty(FluxPotential_solutions{i,j}(1))
                    FPsol(i,j) = FluxPotential_solutions{i,j}(1);
                else
                    error('check');
                end
            end
        end
    end
    %% get contributing rxns
    % we define the Influential Reaction Set (IRS) as the top flux allowance
    % contributors that have a accumulative flux allowance contribution of 95%
    % for each condition, we consider three metrics of the natural distance of
    % all IRS reactions: mean, median, maximum 
    a = find(strcmp(rxnLabel,targetRxn));
    contrRxns = {};
    allowanceTbl = {};
    for j = 1:size(FPsol,2)
        res = FPsol{a,j};
        mytbl = struct();
        mytbl.labels = res.labels;
        mytbl.fluxAllowance = res.full .* res.weight;
        mytbl.UnwDist = res.UnwDist';
        mytbl.flux = res.full;
        mytbl = struct2table(mytbl);
        mytbl = sortrows(mytbl,2,'descend');
        IRS = mytbl.UnwDist(1);
        contrRxns{j,1} = mytbl.labels(1);
        allowanceTbl{j,1} = mytbl(1,:);        
        accFA = mytbl.fluxAllowance(1);
        for k = 2:size(mytbl,1)
           accFA = accFA + mytbl.fluxAllowance(k);
           if accFA < 0.95
               IRS = [IRS,mytbl.UnwDist(k)];
               contrRxns{j,1} = [contrRxns{j,1},mytbl.labels(k)];
               allowanceTbl(j,1) = {[allowanceTbl{j,1};mytbl(k,:)]};
           else
               break;
           end
        end
    end
    rxnSetu = {};
    rxnSeti = contrRxns{1};
    for i = 1:(length(contrRxns)-1)
        rxnSetu = union(rxnSetu,contrRxns{i});
        rxnSeti = intersect(rxnSeti,contrRxns{i});
    end
    % we use the intersecting set of all conditions 
    contributingRxnSets_i{kk} = rxnSeti;
    contributingRxnSets_u{kk} = rxnSetu;
    allowanceTblSets{kk} = allowanceTbl;
end
% %% calculate the coexpression 
% load('output/normalizedLevels_partialExcluded.mat');
% contributingRxnSets = contributingRxnSets_u;
% i = 2;
% targetrxn = targetRxns{i};
% contrRxns = contributingRxnSets{i};
% contrRxns = regexprep(contrRxns,'_f$','');
% contrRxns = regexprep(contrRxns,'_r$','');
% expMat = normalizedLevel_pro_perPro(ismember(valid_rxns_pro_perPro, contrRxns),:);
% corMat = corr(expMat');
% %% make the clustergram and look at co-expressed modules in R!
% histogram(corMat)

%% plot the flux allowance barplots 
i = 2;
% first make the table for the 25 conditions
rxnIDs = contributingRxnSets_u{i};
rxnIDs = regexprep(rxnIDs,'_.$','');
faMat = zeros(length(rxnIDs),25);
for j = 1:size(faMat,2)
    tmp = allowanceTblSets{i}{j};
    [A B] = ismember(contributingRxnSets_u{i},tmp.labels);
    faMat(A,j) = tmp.fluxAllowance(B(A));
end
rxnInfo = rxnIDs;
rxnInfo(:,2) = printRxnFormula_XL(model, rxnIDs,false);
[A B] = ismember(rxnIDs, model.rxns);
rxnInfo(:,3) = model.subSystems(B(A));
for j = 1:length(rxnIDs)
    rxnInfo(j,4) = rxnInfo{j,3}(1);
end
rxnInfo(:,4) = regexprep(rxnInfo(:,4),'^sce.....  ','');
rxnInfo(cellfun(@isempty,rxnInfo(:,4)),4) = {'Not assigned'};
% we manually assign the subsystem for some not assigned ones
rxnInfo(ismember(rxnInfo(:,1),{'r_0454','r_0770','r_0773'}),4) = {'Other energy metabolism'};
rxnInfo(ismember(rxnInfo(:,1),{'r_4214','r_4171'}),4) = {'Cysteine and methionine metabolism'};

% reshape the matrix to make the subsystem-contribution matrix
% subsys = unique(rxnInfo(:,4));
% we only keep a few 
subsys = {'Cysteine and methionine metabolism','Gluconeogenesis','Oxidative phosphorylation','Other energy metabolism','Citrate cycle (TCA cycle)','Others'};
faMat_subsys = zeros(length(subsys),25);
for j = 1:size(faMat_subsys,2)
    for k = 1:length(subsys)
        faMat_subsys(k,j) = sum(faMat(strcmp(subsys{k},rxnInfo(:,4)),j));
    end
end
% reorder
tmp = sum(faMat_subsys,2);
[A B] = sort(tmp,'descend');
faMat_subsys = faMat_subsys(B,:);
subsys = subsys(B);
% assign residuals to others
faMat_subsys(strcmp(subsys,'Others'),:) = 1-sum(faMat_subsys,1);

bar(categorical(conditions),faMat_subsys','stacked')
% rename Gluconeogenesis
subsys(strcmp(subsys,'Gluconeogenesis')) = {'Gluconeogenesis/glycolysis'};
legend(subsys)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5, 4];
plt.LineWidth = 1;
plt.FontSize = 12;
plt.FontName = 'Arial';
ylabel('Flux allowance',  'FontSize',15)
plt.Interpreter = 'none';
plt.LegendLoc = 'eastoutside';
plt.TickLength = [0.01; 0.02];
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.ShowBox = 'off';
plt.LegendLoc = 'eastoutside';
plt.FontSize = 7;
plt.export(['figures/FPA_GST_mechanism_',targetRxns{i},'.pdf']);

%% plot the flux distribution barplots 
% first make the table for the 25 conditions
rxnIDs = contributingRxnSets_u{i};
rxnIDs = regexprep(rxnIDs,'_.$','');
fdMat = zeros(length(rxnIDs),25);
for j = 1:size(fdMat,2)
    tmp = allowanceTblSets{i}{j};
    [A B] = ismember(contributingRxnSets_u{i},tmp.labels);
    fdMat(A,j) = tmp.flux(B(A));
end
rxnInfo = rxnIDs;
rxnInfo(:,2) = printRxnFormula_XL(model, rxnIDs,false);
[A B] = ismember(rxnIDs, model.rxns);
rxnInfo(:,3) = model.subSystems(B(A));
for j = 1:length(rxnIDs)
    rxnInfo(j,4) = rxnInfo{j,3}(1);
end
rxnInfo(:,4) = regexprep(rxnInfo(:,4),'^sce.....  ','');
rxnInfo(cellfun(@isempty,rxnInfo(:,4)),4) = {'Not assigned'};
% we manually assign the subsystem for some not assigned ones
rxnInfo(ismember(rxnInfo(:,1),{'r_0454','r_0770','r_0773'}),4) = {'Other energy metabolism'};
rxnInfo(ismember(rxnInfo(:,1),{'r_4214','r_4171'}),4) = {'Cysteine and methionine metabolism'};

% reshape the matrix to make the subsystem-contribution matrix
% subsys = unique(rxnInfo(:,4));
% we only keep a few 
subsys = {'Cysteine and methionine metabolism','Gluconeogenesis','Oxidative phosphorylation','Other energy metabolism','Citrate cycle (TCA cycle)'};
fdMat_subsys = zeros(length(subsys),25);
for j = 1:size(fdMat_subsys,2)
    for k = 1:length(subsys)
        fdMat_subsys(k,j) = sum(fdMat(strcmp(subsys{k},rxnInfo(:,4)),j));
    end
end
% reorder
tmp = sum(fdMat_subsys,2);
[A B] = sort(tmp,'descend');
fdMat_subsys = fdMat_subsys(B,:);
subsys = subsys(B);

bar(categorical(conditions),fdMat_subsys','stacked')
% rename Gluconeogenesis
subsys(strcmp(subsys,'Gluconeogenesis')) = {'Gluconeogenesis/glycolysis'};
legend(subsys)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [10, 8];
plt.LineWidth = 1;
plt.FontSize = 12;
plt.FontName = 'Arial';
ylabel('Total flux',  'FontSize',15)
plt.Interpreter = 'none';
plt.LegendLoc = 'eastoutside';
plt.TickLength = [0.01; 0.02];
plt.export(['figures/FPA_GST_mechanism_fluxDist_',targetRxns{i},'.tiff']);
%% inspect indiviudal rxns 
% sort the rxn faMat
tmp = sum(faMat,2);
[A B] = sort(tmp,'descend');
faMat = faMat(B,:);
rxnInfo = rxnInfo(B,:);
rxnIDs = rxnIDs(B);
%% load expression data 
load('output/normalizedLevels_partialExcluded.mat');
labels2 = regexprep(conditions,'_','-');

%% plot the expression vs. flux
rxnID = 'r_0892';
printGPRForRxns(model,rxnID);
if any(strcmp(rxnID,valid_rxns_pro_perPro)) 
    myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID),:);
    myFluxLevel = abs(fluxMat_normalized(strcmp(rxnLabel,targetRxn),:));
    figure(7)
    plot(myLevel,myFluxLevel,'.','MarkerSize',10,'Color','k')
    text(myLevel,myFluxLevel,labels2,'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',7)
    xlabel(['Relative expression of ',rxnID])
    ylabel(['Relative flux of ',targetRxn])
    [corrR,corrP] = corr(myFluxLevel',myLevel')
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [3, 2.5];
    plt.LineWidth = 1;
    plt.FontSize = 7;
    plt.XTick = -1:0.2:1;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.Interpreter = 'none';
    plt.ShowBox = 'off';
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.TickDir = 'out';
    plt.export(['figures/FPA_GST_mechanism_expression_',rxnID,'_flux_',targetRxn,'.pdf']);
else
    fprintf('not measured!\n')
end
%% target expression vs flux 
figure;
myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,rxnID),:);
myFlux = abs(fluxMat_normalized(strcmp(rxnLabel,targetRxn),:));
fit = fitlm(myLevel,myFlux);
h = plot(fit);
xlim1 = xlim;
ylim1 = ylim;
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','.')
set(h(1),'MarkerSize',15)
set(h(2), {'color'},{'#808080'}) 
set(h(2), {'LineStyle'},{'--'}) 
set(h(3), {'visible'},{'off'}) 
set(h(4), {'visible'},{'off'}) 
xlabel('Relative Protein Expression');
ylabel('Relative Flux (absolute value)');
text(0.6,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Title = rxnID;
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.XLim = xlim1;
plt.YLim = ylim1;
plt.export(['figures/FPA_GST_mechanism_exp2flux_',rxnID,'.pdf']);
%% plot rFP vs flux
targetRxn = 'r_0485';
dist = 18.5;

load('output/Titration_relativeExp_wtdDist_expDecay_FineGrained.mat');
dorders = n2;
nn = find(n2==dist); 
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

myFlux = abs(fluxMat_normalized(strcmp(rxnLabel,targetRxn),:));
myLevel = relFP(strcmp(valid_rxns,targetRxn),:);

figure;
fit = fitlm(myLevel,myFlux);
h = plot(fit);
xlim1 = xlim;
ylim1 = ylim;
lgd = legend();
set(lgd,'visible','off')
set(h(1), {'color'},{'k'}) 
set(h(1),'Marker','.')
set(h(1),'MarkerSize',15)
set(h(2), {'color'},{'#808080'}) 
set(h(2), {'LineStyle'},{'--'}) 
set(h(3), {'visible'},{'off'}) 
set(h(4), {'visible'},{'off'}) 
xlabel(['rFP of ',targetRxn]);
ylabel(['Relative Flux (absolute value) of ',targetRxn]);
text(0.6,max(myFlux)*0.93,['r = ',num2str(corr(myLevel',myFlux'),2)],'FontSize',15)
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.LineWidth = 1;
plt.FontSize = 7;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.Title = rxnID;
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.XLim = xlim1;
plt.YLim = ylim1;
plt.export(['figures/FPA_GST_flux_rFP_',targetRxn,'_d_',num2str(dist),'.pdf']);
%% finally, robustness analysis of ATP 
model = loadYeatModel();
% the following nutrients need to be set manually
% to start with basic FPA, we allow unlimited exchange
% % phosphate exchange
% model.lb(strcmp(model.rxns,'r_2005')) = r_1244;
% % glucose exchange
% model.lb(strcmp(model.rxns,'r_1714')) = r_1166;
% % ammonium exchange 
% model.lb(strcmp(model.rxns,'r_1654')) = r_1115;
% % uracil 
% model.lb(strcmp(model.rxns,'r_2090')) = r_1272;
% % leucine
% model.lb(strcmp(model.rxns,'r_1899')) = r_1211;
model_ori = model;
fluxMat_raw = fluxMat;
dwTbl = readtable('./../input/YeastJoshua/originalDataTbl/chemostatInfo.xlsx');%gDW/ml cell
factor = repmat(dwTbl.gDCW_mL',size(fluxMat,1),1);
fluxMat_raw = fluxMat_raw * 1000 ./ factor; %mmoles/hr/gDW

rbMat = nan(25,102);
for i = 1:25
    model = model_ori;
    % phosphate exchange
    model.lb(strcmp(model.rxns,'r_2005')) = -fluxMat_raw(strcmp(rxnLabel,'r_1244'),i);
    % glucose exchange
    model.lb(strcmp(model.rxns,'r_1714')) = -fluxMat_raw(strcmp(rxnLabel,'r_1166'),i);
    % ammonium exchange 
    model.lb(strcmp(model.rxns,'r_1654')) = -fluxMat_raw(strcmp(rxnLabel,'r_1115'),i);
    % uracil 
    model.lb(strcmp(model.rxns,'r_2090')) = -fluxMat_raw(strcmp(rxnLabel,'r_1272'),i);
    % leucine
    model.lb(strcmp(model.rxns,'r_1899')) = -fluxMat_raw(strcmp(rxnLabel,'r_1211'),i);
    % maintanence 
    model = changeRxnBounds(model,'r_4046',0,'l'); % maintance 
    model = changeRxnBounds(model,'r_4046',1000,'u'); % maintance 
    % robustness analysis
    model = changeObjective(model,'r_4046');
    maxATP = optimizeCbModel(model);
    maxATP = maxATP.obj;
    model = changeObjective(model,'r_0485');
    for j = 1:101
        model = changeRxnBounds(model,'r_4046',maxATP.*(j-1)/100,'l'); 
        model = changeRxnBounds(model,'r_4046',maxATP.*(j-1)/100,'u'); 
        tmp = optimizeCbModel(model);
        rbMat(i,j) = tmp.obj;
    end
    model = changeRxnBounds(model,'r_4046',0,'l'); % maintance 
    model = changeRxnBounds(model,'r_4046',1000,'u'); % maintance 
    tmp = optimizeCbModel(model);
    rbMat(i,102) = tmp.obj; %no constraint control
end
%% plot
figure
hold on
for i = 1:25
    plot(0:100,rbMat(i,1:101),'-k');
end
hold off
xlabel('ATP maintenance flux (% maximum)');
ylabel('Maximum glutathione synthesis flux (mmoles/gDW/h)');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2, 1.75];
plt.LineWidth = 0.75;
plt.FontSize = 7;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.export(['figures/FPA_GST_glutathione synthesis_flux_robustness.pdf']);

% figure
% hold on
% for i = 1:25
%     plot(0:100,rbMat(i,1:101)./rbMat(i,102),'.-k');
% end
% hold off
