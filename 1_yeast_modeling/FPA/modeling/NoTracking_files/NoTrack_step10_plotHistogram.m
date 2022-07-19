myTbl = readtable('output/coexpression_PCC.csv');
fdr_r = mafdr(myTbl.coexpressionPval,'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
fprintf('%d rxns give significant positive correlation by pearson\n',sum(myTbl.coexpressionPCC(fdr_r<0.05)>0));

figure(1)
hold on
histogram(myTbl.coexpressionPCC,'FaceColor','#0072BD','BinEdges',-1:0.2:1)
xlim([-1,1]);
xlabel('Correlation coefficient');
ylabel('Number of reactions');
histogram(myTbl.coexpressionPCC(fdr_r<0.05 & myTbl.coexpressionPCC >0),'FaceColor','#D95319','BinEdges',-1:0.2:1)
histogram(myTbl.PCC_by_target_expression,'BinEdges',-1:0.2:1)
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
plt.export('figures/coexpression_flux_correlation.pdf');
