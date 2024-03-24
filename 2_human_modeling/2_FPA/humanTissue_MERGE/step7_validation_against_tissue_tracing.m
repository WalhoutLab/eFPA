%% compare with Hui et al 2020 Cell Metab
%% load data 
addpath ../scripts/

% let's start with chow diet 
% error bar ignored - just use mean
TCAcontri_fasted = readtable('input\Hui_et_al_Cell_Metab_2020_dataset\mmc6_formated.xlsx','Sheet','Fasted mice on CD','ReadRowNames',1);
TCAcontri_fed = readtable('input\Hui_et_al_Cell_Metab_2020_dataset\mmc6_formated','Sheet','Fed mice on CD','ReadRowNames',1);
TCAcontri_fasted_KD = readtable('input\Hui_et_al_Cell_Metab_2020_dataset\mmc6_formated.xlsx','Sheet','Fasted mice on KD','ReadRowNames',1);
TCAcontri_fed_KD = readtable('input\Hui_et_al_Cell_Metab_2020_dataset\mmc6_formated','Sheet','Fed mice on KD','ReadRowNames',1);

% averge 
condAve = array2table((TCAcontri_fasted{TCAcontri_fasted_KD.Properties.RowNames,TCAcontri_fasted_KD.Properties.VariableNames} +...
    TCAcontri_fed{TCAcontri_fasted_KD.Properties.RowNames,TCAcontri_fasted_KD.Properties.VariableNames} + ...
    TCAcontri_fasted_KD{TCAcontri_fasted_KD.Properties.RowNames,TCAcontri_fasted_KD.Properties.VariableNames} +...
    TCAcontri_fed_KD{TCAcontri_fasted_KD.Properties.RowNames,TCAcontri_fasted_KD.Properties.VariableNames}) ./ 4);
condAve.Properties.RowNames = TCAcontri_fasted_KD.Properties.RowNames;
condAve.Properties.VariableNames = TCAcontri_fasted_KD.Properties.VariableNames;


% read the FPA results
setEnvForAnalysis
% dataset
load output/FPA_rxn_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
% remove nan
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels(rmInd) = [];
relFP_wtd = relFP;

% normalize to tissue enrichment score
relFP_wtd = relFP_wtd - median(relFP_wtd,2);


% control - no network
load output/FPA_rxn_protein_TS_common_originalFPA_originalDist_order100_naiveNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_noNetwork = [cellfun(@(x) [x,'_f'],targetRxns,'UniformOutput',false);
              cellfun(@(x) [x,'_r'],targetRxns,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_noNetwork(rmInd) = [];
relFP_wtd_noNetwork = relFP;
% normalize to tissue enrichment score
relFP_wtd_noNetwork = relFP_wtd_noNetwork - median(relFP_wtd_noNetwork,2);

%% may be better aligned with simply metabolite degradation potential of a tissue 
% load comsumption potential
% define reaction sets
% the exchange reactions with environment 
excRxns = model.rxns(findExcRxns_XL(model));
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
excRxns = intersect(excRxns,EXinvolvedRxns);

% the transporter (with env) reactions
allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
TSP = [];
for i = 1:length(allCmp_iHumanName)
    myMet_e = {[allCmp_iHumanName{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,allCmp_iHumanName(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    % we define the transporter as the reactions that contain the same
    % metabolite in [e] and another compartment (cellular) in the same
    % reaction
    candidate = intersect(myRxns_non_e,myRxns_e);
    % check if is on diff side of the reaction
    if ~isempty(candidate)
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                TSP = union(TSP,candidate(j));
            end
        end
    end
end
% some special transporter will be missed, we add back 
envTspRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
envTspRxns = intersect(envTspRxns,EXinvolvedRxns);
tspRxns = union(TSP, envTspRxns);

% the internal transporters 
% internal transporters are hard to define, since a reaction can span two
% compartments internally but it is not a real transporter. So, we use the
% subsys annotation as a compromise
intTspRxns = setdiff(model.rxns(ismember(model.subSystems,{'Transport reactions'})),tspRxns);

% internal regular reactions
intRxns = setdiff(model.rxns, [excRxns; tspRxns; intTspRxns]);

% regular met-analysis target (transporters)
targetRxns = tspRxns;

load('input/allCmp_iHumanName.mat');
allCmp_iHumanName = unique(allCmp_iHumanName);
targetRxns_DM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);

allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
targetRxns_allMetDM = cellfun(@(x) ['NewMet_',x],allCmp_iHumanName,'UniformOutput',false);

load output/FPA_allMetDemand_protein_TS_common_newFPA_weightedDist_order6_tissueNetwork.mat
relFP = [relFP_f;relFP_r];
rowlabels_allMetDM = [cellfun(@(x) [x,'_f'],targetRxns_allMetDM,'UniformOutput',false);
                cellfun(@(x) [x,'_r'],targetRxns_allMetDM,'UniformOutput',false)];
rmInd = all(isnan(relFP),2);
relFP(rmInd,:) = [];
rowlabels_allMetDM(rmInd) = [];
relFP_wtd_allMetDM = relFP;

% normalize to tissue enrichment score
relFP_wtd_allMetDM = relFP_wtd_allMetDM - median(relFP_wtd_allMetDM,2);


%% glucose usage 

% not very interesting because the entire pathway is very well coexpressed
% and the flux can be predicted just by the expression of any gene 

glucose2tca = {'HMR_4394_f'};%HMR_4394_f,'HMR_4381_r','HMR_4379_f','HMR_4375_r','HMR_4373_r','HMR_4368_f','HMR_4365_r','HMR_4363_f','HMR_4358_f','HMR_4137_f'};

% we use the hxk reaction as a prediction readout 
hxk_FPA = mean(relFP_wtd(ismember(rowlabels,glucose2tca),:),1); % 

% align the tissues 
tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartVentricle', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial
tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};

% plot data 
[A B] = ismember(tissue_pred, conditions);
FPA_data = hxk_FPA(B(A));
[A B] = ismember(tissue_data, condAve.Properties.RowNames);
tracing_data = condAve.Glucose(B(A))';

% Create the scatter plot
figure;
scatter(FPA_data, tracing_data, 'filled');
xlabel('delta rFP of glycolysis-PDH pathway reactions');
ylabel('Fractional contribution of circulation glucose to tissue TCA cycle');

% Add labels
for i = 1:length(FPA_data)
    text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

[r p] =corr(FPA_data', tracing_data','type','Pearson');
title(['PCC=',num2str(r),', P=',num2str(p)])


plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2*1.5, 1.75*1.5];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';

plt.export(['figures/validation_with_public_tracing_data_glucose_usage.pdf']);


%% glucose usage - controls - no integration

% looks fed state gives better alignment - since we cannot predict
% different states, we just arbiturily go with one for presentation

glucose2tca = {'HMR_4394_f'};%,'HMR_4381_r','HMR_4379_f','HMR_4375_r','HMR_4373_r','HMR_4368_f','HMR_4365_r','HMR_4363_f','HMR_4358_f','HMR_4137_f'};

% we use the hxk reaction as a prediction readout 
hxk_FPA = mean(relFP_wtd_noNetwork(ismember(rowlabels_noNetwork,glucose2tca),:),1); % 

% align the tissues 
tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartVentricle', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial
tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};

% plot data 
[A B] = ismember(tissue_pred, conditions);
FPA_data = hxk_FPA(B(A));
[A B] = ismember(tissue_data, condAve.Properties.RowNames);
tracing_data = condAve.Glucose(B(A))';

% Create the scatter plot
figure;
scatter(FPA_data, tracing_data, 'filled');
xlabel('delta rFP of glycolysis-PDH pathway reactions');
ylabel('Fractional contribution of circulation glucose to tissue TCA cycle');

% Add labels
for i = 1:length(FPA_data)
    text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

[r p] =corr(FPA_data', tracing_data','type','Pearson');
title(['PCC=',num2str(r),', P=',num2str(p)])


plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2*1.5, 1.75*1.5];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';

plt.export(['figures/validation_with_public_tracing_data_glucose_usage_no_integration_control.pdf']);

%% peroxi beta oxi vignette 

% we use the total FA contribution to TCA as a proxy of the peroxi beta oxi
% flux 

% we use the hxk reaction as a prediction readout FAOXC160_fz
hxk_FPA = mean(relFP_wtd(ismember(rowlabels,{'HMR_3078_f'}),:),1); % FAOXC160_f
% align the tissues 
tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartAtrial', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial HeartVentricle
tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};

% plot data 
[A B] = ismember(tissue_pred, conditions);
FPA_data = hxk_FPA(B(A));
[A B] = ismember(tissue_data, condAve.Properties.RowNames);
tracing_data = sum(condAve{B(A), {'OleicAcid','LinoleicAcid','PalmiticAcid'}},2)';

% Create the scatter plot
figure;
scatter(FPA_data, tracing_data, 'filled');
xlabel('delta rFP of peroxi-beta oxidation reaction (HMR\_3078)');
ylabel('Total fractional contribution of circulation FAs to tissue TCA cycle');

% Add labels
for i = 1:length(FPA_data)
    text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

[r p] =corr(FPA_data', tracing_data','type','Pearson');
title(['PCC=',num2str(r),', P=',num2str(p)])

plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2*1.5, 1.75*1.5];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.export(['figures/validation_with_public_tracing_data_FA_peroxi_oxi_usage.pdf']);

% we use the total FA contribution to TCA as a proxy of the peroxi beta oxi
% flux 

% we use the hxk reaction as a prediction readout FAOXC160_fz
hxk_FPA = mean(relFP_wtd_noNetwork(ismember(rowlabels_noNetwork,{'HMR_3078_f'}),:),1); % FAOXC160_f
% align the tissues 
tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartAtrial', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial HeartVentricle
tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};

% plot data 
[A B] = ismember(tissue_pred, conditions);
FPA_data = hxk_FPA(B(A));
[A B] = ismember(tissue_data, condAve.Properties.RowNames);
tracing_data = sum(condAve{B(A), {'OleicAcid','LinoleicAcid','PalmiticAcid'}},2)';

% Create the scatter plot
figure;
scatter(FPA_data, tracing_data, 'filled');
xlabel('delta rFP (no-integration) of peroxi-beta oxidation reaction (HMR\_3078)');
ylabel('Total fractional contribution of circulation FAs to tissue TCA cycle');

% Add labels
for i = 1:length(FPA_data)
    text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

[r p] =corr(FPA_data', tracing_data','type','Pearson');
title(['PCC=',num2str(r),', P=',num2str(p)])

plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [2*1.5, 1.75*1.5];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XTick = -1:0.2:1;
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.ShowBox = 'off';
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickDir = 'out';
plt.Interpreter = 'None';
plt.export(['figures/validation_with_public_tracing_data_FA_peroxi_oxi_usage_no_integration_control.pdf']);

%% lactate usage 
% lactate usage is not well predicted but likely because exp. contribution
% of lactate is always significant in most tissues
% 
% lact2tca = {'HMR_4388_r','HMR_4137_f'};
% 
% we use the hxk reaction as a prediction readout 
% hxk_FPA = mean(relFP_wtd(ismember(rowlabels,lact2tca),:),1); % 
% 
% align the tissues 
% tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartVentricle', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial
% tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};
% 
% plot data 
% [A B] = ismember(tissue_pred, conditions);
% FPA_data = hxk_FPA(B(A));
% [A B] = ismember(tissue_data, TCAcontri_fed.Properties.RowNames);
% tracing_data = TCAcontri_fed.Lactate(B(A))';
% 
% Create the scatter plot
% figure;
% scatter(FPA_data, tracing_data, 'filled');
% xlabel('delta rFP of LDH-PDH reactions');
% ylabel('Fractional contribution of circulation glucose to tissue TCA cycle');
% 
% Add labels
% for i = 1:length(FPA_data)
%     text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end
% 
% [r p] =corr(FPA_data', tracing_data','type','Pearson');
% title(['PCC=',num2str(r),', P=',num2str(p)])
% 
% 
% % Alanine usage 
% also not good; maybe only if the pathway can be unambiguously defined it
% could work
% 
% lact2tca = {'HMR_4109_f','HMR_4137_f'};
% 
% we use the hxk reaction as a prediction readout 
% hxk_FPA = mean(relFP_wtd(ismember(rowlabels,lact2tca),:),1); % 
% 
% align the tissues 
% tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartVentricle', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial
% tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};
% 
% plot data 
% [A B] = ismember(tissue_pred, conditions);
% FPA_data = hxk_FPA(B(A));
% [A B] = ismember(tissue_data, TCAcontri_fed.Properties.RowNames);
% tracing_data = TCAcontri_fed.Alanine(B(A))';
% 
% Create the scatter plot
% figure;
% scatter(FPA_data, tracing_data, 'filled');
% xlabel('delta rFP of LDH-PDH reactions');
% ylabel('Fractional contribution of circulation glucose to tissue TCA cycle');
% 
% Add labels
% for i = 1:length(FPA_data)
%     text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end
% 
% [r p] =corr(FPA_data', tracing_data','type','Pearson');
% title(['PCC=',num2str(r),', P=',num2str(p)])
% 
% % glutamine usage - not very well - but the route is hard to define 
% 
% looks fed state gives better alignment - since we cannot predict
% different states, we just arbiturily go with one for presentation
% 
% we use the hxk reaction as a prediction readout 
% hxk_FPA = relFP_wtd(strcmp(rowlabels,'HMR_4197_f'),:); % AADAT reaction
% align the tissues 
% tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartVentricle', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial
% tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};
% 
% plot data 
% [A B] = ismember(tissue_pred, conditions);
% FPA_data = hxk_FPA(B(A));
% [A B] = ismember(tissue_data, TCAcontri_fed.Properties.RowNames);
% tracing_data = TCAcontri_fed.Glutamine(B(A))';
% 
% Create the scatter plot
% figure;
% scatter(FPA_data, tracing_data, 'filled');
% xlabel('delta rFP of AADAT reaction');
% ylabel('Fractional contribution of circulation glucose to tissue TCA cycle');
% 
% Add labels
% for i = 1:length(FPA_data)
%     text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
% end


%% BCAA usage - not very well - but the route is hard to define 

% looks fed state gives better alignment - since we cannot predict
% different states, we just arbiturily go with one for presentation

% we use the hxk reaction as a prediction readout 
hxk_FPA = mean(relFP_wtd_noNetwork(ismember(rowlabels_noNetwork,{'HMR_3744_f','HMR_3765','HMR_3778'}),:),1); % BCAT2 reaction in valine deg 
% align the tissues 
tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartVentricle', 'SmallIntestine'};
% try on BrainCerebellum HeartAtrial
tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};

% plot data 
[A B] = ismember(tissue_pred, conditions);
FPA_data = hxk_FPA(B(A));
[A B] = ismember(tissue_data, TCAcontri_fasted.Properties.RowNames);
tracing_data = mean([TCAcontri_fasted.BCAAs(B(A)),TCAcontri_fed.BCAAs(B(A))],2);

% Create the scatter plot
figure;
scatter(FPA_data, tracing_data, 'filled');
xlabel('delta rFP of BCAT2 reaction');
ylabel('Fractional contribution of circulation glucose to tissue TCA cycle');

% Add labels
for i = 1:length(FPA_data)
    text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

[r p] =corr(FPA_data', tracing_data','type','Pearson');
title(['PCC=',num2str(r),', P=',num2str(p)])


%% directly use metabolite comsumption potential 
% looks fed state gives better alignment - since we cannot predict
% different states, we just arbiturily go with one for presentation

metName1_list = {'glucose','L-lactate','glutamine','alanine','(R)-3-hydroxybutanoate','palmitate','oleate',...
    'linoleate','serine','glycine','citrate','acetate',{'leucine','isoleucine','valine'},{'methionine','lysine','phenylalanine','tryptophan','threonine','histidine'}};
metName2_list = TCAcontri_fed.Properties.VariableNames;

for zz = 1:length(metName1_list) 
    metName1 = metName1_list{zz}; % {'leucine','isoleucine','valine'};
    metName2 = metName2_list{zz}; % 'BCAAs';
    
    % we use the hxk reaction as a prediction readout 
    met_cons_FPA = mean(relFP_wtd_allMetDM(ismember(rowlabels_allMetDM,strcat('NewMet_',metName1,'_r')),:),1);
    
    % align the tissues 
    tissue_pred = {'BrainCortex', 'MuscleSkeletal', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'HeartVentricle', 'SmallIntestine'};
    % try on BrainCerebellum HeartAtrial
    tissue_data = {'Brain', 'Muscle', 'Liver', 'Spleen', 'Pancreas', 'Lung', 'Heart', 'Small intestine'};
    
    % plot data 
    [A B] = ismember(tissue_pred, conditions);
    FPA_data = met_cons_FPA(B(A));
    [A B] = ismember(tissue_data, TCAcontri_fed.Properties.RowNames);
    tracing_data = TCAcontri_fed.(metName2)(B(A))';
    
    % Create the scatter plot
    figure;
    scatter(FPA_data, tracing_data, 'filled');
    xlabel(['delta consumption potential of ',metName1]);
    ylabel(['Fractional contribution of circulation ',metName2,' to tissue TCA cycle']);
    
    % Add labels
    for i = 1:length(FPA_data)
        text(FPA_data(i), tracing_data(i), tissue_data{i}, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    end
    
    [r p] =corr(FPA_data', tracing_data');
    title(['PCC=',num2str(r),', P=',num2str(p)])
end
