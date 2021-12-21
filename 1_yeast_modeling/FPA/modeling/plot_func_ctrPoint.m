function wd = plot_func_ctrPoint(name,ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model) 
wd = ['figures/controlingPoints/',name];
mkdir(wd);
delete([wd,'/RxnInfo.txt']);
diary([wd,'/RxnInfo.txt']);
fprintf('controlling point:\n');
printRxnFormula_XL(model,ctrPoint);
fprintf('\ntarget rxns:\n');
printRxnFormula_XL(model,controlledRxns);
diary off
for i = 1:length(controlledRxns)
    myLevel = normalizedLevel_pro_perPro(strcmp(valid_rxns_pro_perPro,ctrPoint),:);
    myFluxLevel = abs(fluxMat_normalized(strcmp(rxnLabel,controlledRxns{i}),:));
    figure(1)
    fit = fitlm(myLevel,myFluxLevel);
    [r p] = corr(myLevel',myFluxLevel');
    h = plot(fit);
    lgd = legend();
    set(lgd,'visible','off')
    set(h(1), {'color'},{'k'}) 
    set(h(1),'Marker','.')
    set(h(1),'MarkerSize',15)
    set(h(2), {'color'},{'#808080'}) 
    set(h(2), {'LineStyle'},{'--'}) 
    set(h(3), {'visible'},{'off'}) 
    set(h(4), {'visible'},{'off'}) 
    xlabel(['Relative expression of ',ctrPoint])
    ylabel(['Relative flux of ',controlledRxns{i}])
    text(0.47,max(myFluxLevel)*0.93,['r = ',num2str(r,2)],'FontSize',7)
    text(0.47,max(myFluxLevel)*0.85,['p < ',num2str(p,2)],'FontSize',7)
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [2, 1.75];
    plt.LineWidth = 1;
    plt.FontSize = 7;
    % plt.XTick = -1:0.2:1;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.ShowBox = 'off';
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.TickDir = 'out';
    plt.Interpreter = 'None';
    plt.export([wd, '/flux_exp_correlation_ctrPoint_',ctrPoint,'_targetRxn_',controlledRxns{i},'.pdf']);
end
%% plot flux vs. flux correlation matrix
for i = 1:length(controlledRxns)
    myLevel = abs(fluxMat_normalized(strcmp(rxnLabel,ctrPoint),:));
    myFluxLevel = abs(fluxMat_normalized(strcmp(rxnLabel,controlledRxns{i}),:));
    figure(1)
    fit = fitlm(myLevel,myFluxLevel);
    [r p] = corr(myLevel',myFluxLevel');
    h = plot(fit);
    lgd = legend();
    set(lgd,'visible','off')
    set(h(1), {'color'},{'k'}) 
    set(h(1),'Marker','.')
    set(h(1),'MarkerSize',15)
    set(h(2), {'color'},{'#808080'}) 
    set(h(2), {'LineStyle'},{'--'}) 
    set(h(3), {'visible'},{'off'}) 
    set(h(4), {'visible'},{'off'}) 
    xlabel(['Relative flux of control point (',ctrPoint,')'])
    ylabel(['Relative flux of target (',controlledRxns{i},')'])
    text(max(myLevel)*0.5,max(myFluxLevel)*0.93,['r = ',num2str(r,2)],'FontSize',7)
    text(max(myLevel)*0.5,max(myFluxLevel)*0.85,['p < ',num2str(p,2)],'FontSize',7)
    plt = Plot(); % create a Plot object and grab the current figure
    plt.BoxDim = [2, 1.75];
    plt.LineWidth = 1;
    plt.FontSize = 7;
    plt.LegendLoc = 'NorthWest';
    plt.FontName = 'Arial';
    plt.ShowBox = 'off';
    plt.XMinorTick = 'off';
    plt.YMinorTick = 'off';
    plt.TickDir = 'out';
    plt.Interpreter = 'None';
    plt.export([wd, '/flux_flux_correlation_ctrPoint_',ctrPoint,'_targetRxn_',controlledRxns{i},'.pdf']);
end
