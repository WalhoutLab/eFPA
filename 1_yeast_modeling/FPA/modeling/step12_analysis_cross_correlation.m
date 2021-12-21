%% rationale
% systematically check the cross-correlation between the flux of ROI and
% any measured expression (at reaction level). This may indicate controling
% point in a pathway, or distal relations between two reactions. 
% see if this search will get more flux explained (>44%) or less, or bring
% something new.
%% load the model
addpath('./../scripts/')
model = loadYeatModel();
%% load the expression files and other optional inputs
% load expression data 
proTbl = readtable('./../input/YeastJoshua/originalDataTbl/proteinTbl.xlsx');% this is the log2(FC_reference)
conditions = proTbl.Properties.VariableNames(2:end);
load('output/normalizedLevels_partialExcluded.mat')
normalizedLevel = normalizedLevel_pro_perPro;
valid_rxns = valid_rxns_pro_perPro ;
% flux table 
fluxTbl = readtable('./../input/YeastJoshua/originalDataTbl/fluxTbl.xlsx');
% make the matched matrix
fluxMat = [];
for i = 1: length(conditions)
    fluxMat(:,i) = (fluxTbl.([conditions{i},'_FVA_MAX']) + fluxTbl.([conditions{i},'_FVA_MIN']))/2;
end
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat = fluxMat ./ repmat(GRrate.DR_Actual',size(fluxMat,1),1);
fluxMat_normalized = fluxMat;
rxnLabel = fluxTbl.Model_Reaction_ID;
%% Computing the correlation between any expression/reaction pair 
% we construct the correlation metric matrix: flux (rows) -> expression
% (columns)
r = zeros(length(rxnLabel),length(valid_rxns));
p_r = ones(length(rxnLabel),length(valid_rxns));
for j = 1:length(valid_rxns)
    myExp = normalizedLevel(j,:)';
    for i = 1:length(rxnLabel)
        myFlux = abs(fluxMat(i,:))';
        [a1,a2] = corr(myExp,myFlux,'type','Pearson');
        r(i,j) = a1;
        p_r(i,j) = a2;
    end
end
% calculate the FDR matrix 
FDRr = ones(length(rxnLabel),length(valid_rxns));
FDRc = ones(length(rxnLabel),length(valid_rxns));
% we did the multiple test in two-dimension, so we had to correct for both 
for i = 1:length(rxnLabel) % for each reaction, we correct for multiple testing against screening all expression vectors
    FDRr(i,:) = mafdr(p_r(i,:),'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
end
for i = 1:length(valid_rxns) % for each reaction, we correct for multiple testing against screening all expression vectors
    FDRc(:,i) = mafdr(p_r(:,i),'BHFDR', true);% BHFDR adjustment is chosen to keep strigency (control FDR instead of estimate FDR (which is used in pFDR(qvalue)))
end
FDR = 1-(1-FDRr).*(1-FDRc);

% histogram(max(r,[],2))
sum(any(r > 0 & FDR < 0.05, 2))
figure;
hold on
histogram(max(r,[],2));
histogram(max(r(any(r > 0 & FDR < 0.05, 2),:),[],2))
hold off 

% save data and perform further visualization & analysis in R
t = array2table(r);
t.Properties.RowNames = rxnLabel;
t.Properties.VariableNames = valid_rxns;
writetable(t,'output/crossCorrelation_rMat.csv','WriteRowNames',1);
t = array2table(FDR);
t.Properties.RowNames = rxnLabel;
t.Properties.VariableNames = valid_rxns;
writetable(t,'output/crossCorrelation_FDRMat.csv','WriteRowNames',1);

%% plot the flux-flux-flux-expression co-correlation matrix for locating potential controling points
myRxns = intersect(rxnLabel,valid_rxns);
flux_flux_corr = corr(abs(fluxMat)'); % the direction of flux mean be meaningless, so we use absolute level 
[A B] = ismember(myRxns,rxnLabel);
flux_flux_corr = flux_flux_corr(B(A), B(A));
flux_exp_corr = r(B(A), :);
[A B] = ismember(myRxns,valid_rxns);
flux_exp_corr = flux_exp_corr(:, B(A));
coCorMat = flux_flux_corr .* flux_exp_corr;
% save data and perform further visualization & analysis in R
t = array2table(coCorMat);
t.Properties.RowNames = myRxns;
t.Properties.VariableNames = myRxns;
writetable(t,'output/coCorrelation_rMat.csv','WriteRowNames',1);

% we make a strigent set that also requires self correlation
selfcorr = diag(flux_exp_corr)';
selfcorr(selfcorr<0) = 0; % only pos self corr is meaningful
coCorMat2 = flux_flux_corr .* flux_exp_corr .* repmat(selfcorr,size(flux_exp_corr,1),1);
% save data and perform further visualization & analysis in R
t = array2table(coCorMat2);
t.Properties.RowNames = myRxns;
t.Properties.VariableNames = myRxns;
writetable(t,'output/coCorrelation_rMat_stringent.csv','WriteRowNames',1);


%% check some possible network features of correlated rxns
%% compare the distance with correlation 
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
distMat_wtd = distMat_min;
labels_ID = regexprep(labels,'_.$','');
a = zeros(length(model.rxns),length(model.rxns));
inds = false(length(model.rxns),length(labels_ID));
for j = 1:length(model.rxns)
    inds(j,:) = strcmp(labels_ID, model.rxns{j});
end
for i = 1:length(model.rxns)
    ind = inds(i,:);
    if any(ind)
        tmp = min(distMat_min(ind,:),[],1);
        for j = 1:length(model.rxns)
            ind = inds(j,:);
            if any(ind)
                a(i,j) = min(tmp(ind));
            end
        end
    end
    i
end

distance_raw = readtable('./../input/YeastJoshua/distanceMatrix.txt','FileType','text','ReadRowNames',true); %we load from the output of the distance calculator. For usage of distance calculator, please refer to the section in Github
labels = distance_raw.Properties.VariableNames;
labels = cellfun(@(x) [x(1:end-1),'_',x(end)],labels,'UniformOutput',false);
distMat_raw = table2array(distance_raw);
distMat_min = zeros(size(distance_raw,1),size(distance_raw,2));
for i = 1:size(distMat_min,1)
    for j = 1:size(distMat_min,2)
        distMat_min(i,j) = min([distMat_raw(i,j),distMat_raw(j,i)]);
    end
end
distMat_normal = distMat_min;
labels_ID = regexprep(labels,'_.$','');
a = zeros(length(model.rxns),length(model.rxns));
inds = false(length(model.rxns),length(labels_ID));
for j = 1:length(model.rxns)
    inds(j,:) = strcmp(labels_ID, model.rxns{j});
end
for i = 1:length(model.rxns)
    ind = inds(i,:);
    if any(ind)
        tmp = min(distMat_min(ind,:),[],1);
        for j = 1:length(model.rxns)
            ind = inds(j,:);
            if any(ind)
                a(i,j) = min(tmp(ind));
            end
        end
    end
    i
end
% map to the r_mat
[A1 B1] = ismember(valid_rxns, model.rxns);
[A2 B2] = ismember(rxnLabel, model.rxns);
dMat_normal = distMat_normal(B2(A2), B1(A1));
dMat_wtd = distMat_wtd(B2(A2), B1(A1));
dMat_normal(isinf(dMat_normal)) = max(max(dMat_normal(~isinf(dMat_normal))));
dMat_wtd(isinf(dMat_wtd)) = max(max(dMat_wtd(~isinf(dMat_wtd))));
dMat_normal(isnan(dMat_normal)) = max(max(dMat_normal(~isnan(dMat_normal))));
dMat_wtd(isnan(dMat_wtd)) = max(max(dMat_wtd(~isnan(dMat_wtd))));

plot(r(:), (dMat_normal(:)),'.','MarkerSize',0.1)
plot(r(:), (dMat_normal(:)),'.','MarkerSize',0.1)

NoTrack_temp_scartPlot(r(~isnan(dMat_wtd)), (dMat_wtd(~isnan(dMat_wtd))))

NoTrack_temp_scartPlot(r(1,~isnan(dMat_wtd(1,:))), (dMat_wtd(1,~isnan(dMat_wtd(1,:)))))

%% check for the distribution of distance for each level of correlation 
figure
hold on 
h1=histogram(dMat_normal,'Normalization','probability');
histogram(dMat_normal(r>0 & FDR <0.05),'Normalization','probability','BinEdges',h1.BinEdges);
legend({'ctr','sig'})
hold off
%%
figure
hold on 
h1=histogram(dMat_wtd,'Normalization','probability');
histogram(dMat_wtd(r>0 & FDR <0.05),'Normalization','probability','BinEdges',h1.BinEdges);
legend({'ctr','sig'})
hold off

% ==> didnt find obvious feature in the network distance (between ROI and
% correlated enzyme) but noted that ~25% of correlated enzyme is not
% connected to ROI (inf distance), this may indicate the correlation is a
% false discovery 




