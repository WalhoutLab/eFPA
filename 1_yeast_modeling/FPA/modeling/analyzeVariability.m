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

load('output/normalizedLevels_partialExcluded.mat');
expMat_normalized = normalizedLevel_pro_perPro(ismember(valid_rxns_pro_perPro,rxnLabel),:);
%% find jumps in flux and expression matrix
AUC_flux = [];
AUC_exp = [];
k0 = 0.66;
isChange_flux = min(fluxMat_normalized,[],2)<k0;
sum(isChange_flux)
isChange_exp = min(expMat_normalized,[],2)<k0;
sum(isChange_exp)

for i = 1:size(fluxMat_normalized)
    sorted = sort(fluxMat_normalized(i,:));
    % scale to 0-1
    sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
    steps = sorted(2:end) - sorted(1:end-1);
    steps_sorted = sort(steps); 
    % reordered curve is
    reordered = [];
    reordered(1) = 0;
    for z = 2:length(sorted)
        reordered(z) = sum(steps_sorted(1:z-1));
    end
    AUC_flux(i) = trapz(1:25,reordered)/12;
end

for i = 1:size(expMat_normalized)
    sorted = sort(expMat_normalized(i,:));
    % scale to 0-1
    sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
    steps = sorted(2:end) - sorted(1:end-1);
    steps_sorted = sort(steps); 
    % reordered curve is
    reordered = [];
    reordered(1) = 0;
    for z = 2:length(sorted)
        reordered(z) = sum(steps_sorted(1:z-1));
    end
    AUC_exp(i) = trapz(1:25,reordered)/12;
end
    
figure;
hold on
histogram(AUC_flux);
histogram(AUC_exp);
hold off
histogram([AUC_exp,AUC_flux]);

AUCk = quantile([AUC_exp,AUC_flux],0.25);
%% visualize seperation and correlation of expression
isjump_flux = false(size(fluxMat_normalized,1),1);
isjump_exp = false(size(expMat_normalized,1),1);
for i = 1:size(expMat_normalized)
    sorted = sort(expMat_normalized(i,:));
    % scale to 0-1
    sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
    steps = sorted(2:end) - sorted(1:end-1);
    steps_sorted = sort(steps); 
    % reordered curve is
    reordered = [];
    reordered(1) = 0;
    for z = 2:length(sorted)
        reordered(z) = sum(steps_sorted(1:z-1));
    end
    AUC = trapz(1:25,reordered)/12;
    if AUC < AUCk
        isjump_exp(i) = true;
    else
        isjump_exp(i) = false;
    end
end
for i = 1:size(fluxMat_normalized)
    sorted = sort(fluxMat_normalized(i,:));
    % scale to 0-1
    sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
    steps = sorted(2:end) - sorted(1:end-1);
    steps_sorted = sort(steps); 
    % reordered curve is
    reordered = [];
    reordered(1) = 0;
    for z = 2:length(sorted)
        reordered(z) = sum(steps_sorted(1:z-1));
    end
    AUC = trapz(1:25,reordered)/12;
    if AUC < AUCk
        isjump_flux(i) = true;
    else
        isjump_flux(i) = false;
    end
end

% plot
jumpping_flux = fluxMat_normalized(isjump_flux&isChange_flux,:);
jumpping_exp = expMat_normalized(isjump_exp&isChange_exp,:);
nojumpping_flux = fluxMat_normalized(~isjump_flux&isChange_flux,:);
nojumpping_exp = expMat_normalized(~isjump_exp&isChange_exp,:);figure;
figure
hold on 
% scale to 0-1
nojumpping_flux = (nojumpping_flux - min(nojumpping_flux,[],2)) ./ repmat((max(nojumpping_flux,[],2) - min(nojumpping_flux,[],2)),1,25);
for i = 1:size(nojumpping_flux,1)
    plot(sort(nojumpping_flux(i,:)),'g');
end
jumpping_flux = (jumpping_flux - min(jumpping_flux,[],2)) ./ repmat((max(jumpping_flux,[],2) - min(jumpping_flux,[],2)),1,25);
for i = 1:size(jumpping_flux,1)
    plot(sort(jumpping_flux(i,:)),'r');
end
hold off
xlabel('rank');
ylabel('scaled flux');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.7, 4.7];
plt.LineWidth = 0.5;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.LegendLoc = 'northwest';
plt.export(['figures/show_flux_variability.tiff']);
% scale to 0-1
figure
hold on
nojumpping_exp = (nojumpping_exp - min(nojumpping_exp,[],2)) ./ repmat((max(nojumpping_exp,[],2) - min(nojumpping_exp,[],2)),1,25);
for i = 1:size(nojumpping_exp,1)
    plot(sort(nojumpping_exp(i,:)),'g');
end
jumpping_exp = (jumpping_exp - min(jumpping_exp,[],2)) ./ repmat((max(jumpping_exp,[],2) - min(jumpping_exp,[],2)),1,25);
for i = 1:size(jumpping_exp,1)
    plot(sort(jumpping_exp(i,:)),'r');
end
hold off
xlabel('rank');
ylabel('scaled expression');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [5.7, 4.7];
plt.LineWidth = 0.5;
plt.FontSize = 15;
plt.FontName = 'Arial';
plt.LegendLoc = 'northwest';
plt.export(['figures/show_expression_variability.tiff']);


%%  seperate out the jummping fluxes (for precision)
