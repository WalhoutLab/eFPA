figure;
hold on
myRxns = sigCorr.rxnID(strcmp(sigCorr.correlated,'Yes'));
[A B] = ismember(myRxns, rxnLabel);
flux_flux_corr_subset = flux_flux_corr(B(A),:);

nRxns = [];
for i = 0:0.01:1
    nRxns = [nRxns,sum(max(abs(flux_flux_corr_subset),[],1) > i)]; % or 0.99 (which gives 88 vs. 161 here)
end

plot(0:0.01:1,nRxns,'-r')
for k = 1:100
    myRxns = sigCorr.rxnID(randperm(length(sigCorr.rxnID),46));
    [A B] = ismember(myRxns, rxnLabel);
    flux_flux_corr_subset = flux_flux_corr(B(A),:);

    nRxns = [];
    for i = 0:0.01:1
        nRxns = [nRxns,sum(max(abs(flux_flux_corr_subset),[],1) > i)]; % or 0.99 (which gives 88 vs. 161 here)
    end

    plot(0:0.01:1,nRxns,'-k')
end
hold off