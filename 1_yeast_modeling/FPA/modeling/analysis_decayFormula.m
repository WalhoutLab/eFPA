%%
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
distMat = distMat_min;
figure
histogram(distMat);
ylabel('Number of reaction pairs');
xlabel('Metabolic distance');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3, 2.4];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/distance_distance_distribution_unweighted.pdf']);
%%
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
figure;
histogram(distMat,'BinEdges',0:1:65);
ylabel('Number of reaction pairs');
xlabel('Metabolic distance');
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [3, 2.4];
plt.LineWidth = 1;
plt.FontSize = 10;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.LegendLoc = 'NorthWest';
plt.FontName = 'Arial';
plt.export(['figures/distance_distance_distribution_weighted.pdf']);

%%
d = 0:0.5:20;
colorList = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','#FF00FF','#FF0000',[1 1 1].*0.5,'k'};
%% order law decay
clear y
n = 0:1:10;
for i = 1:length(n)
    y(i,:) = 100 ./ (1+d).^n(i);
end
figure;
hold on 
for i = 1:length(n)
    plot(d,y(i,:),'.-','Color',colorList{i});
end
hold off
xlabel('Metabolic distance')
ylabel('Weight decay level (original weight%)');
legend(cellfun(@(x) ['order = ',x], strsplit(num2str(n)),'UniformOutput',0))
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.YLim = [0 105];
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.LegendLoc = 'East';
%plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.export(['figures/distance_orderLaw_decay.pdf']);
%% exponential decay - base 2
figure
k = 2;
clear y
n = 0:1:10;
for i = 1:length(n)
    y(i,:) = 100 ./ (1+k.^(d-n(i)));
end
figure(1)
hold on 
for i = 1:length(n)
    plot(d,y(i,:),'-','Color',colorList{i});
end
hold off
xlabel('Metabolic distance')
ylabel('Weight decay level (original weight%)');
legend(cellfun(@(x) ['dist. bound. = ',x], strsplit(num2str(n)),'UniformOutput',0))
title('exponential base = 2')
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.YLim = [0 105];
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.LegendLoc = 'NorthEast';
%plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.export(['figures/distance_orderLaw_decay_base2.pdf']);
%% exponential decay - base e
% clear y
% n = 0:1:10;
% for i = 1:length(n)
%     y(i,:) = 100 ./ (1+exp(d-n(i)));
% end
% figure(1)
% hold on 
% for i = 1:length(n)
%     plot(d,y(i,:),'.-');
% end
% hold off
% xlabel('distance')
% ylabel('dilution level (% of original weight)');
% legend(strsplit(num2str(n)))
% title('decay by base e')
%% exponential decay - base 100
figure
k = 1000;
clear y
n = 0:1:10;
for i = 1:length(n)
    y(i,:) = 100 ./ (1+k.^(d-n(i)));
end
figure(1)
hold on 
for i = 1:length(n)
    plot(d,y(i,:),'-','Color',colorList{i});
end
hold off
xlabel('Metabolic distance')
ylabel('Weight decay level (original weight%)');
legend(cellfun(@(x) ['dist. bound. = ',x], strsplit(num2str(n)),'UniformOutput',0))
title('exponential base = 1000')
plt = Plot(); % create a Plot object and grab the current figure
plt.BoxDim = [1.95, 1.6125];
plt.LineWidth = 1;
plt.YLim = [0 105];
plt.FontSize = 7;
plt.FontName = 'Arial';
plt.LegendLoc = 'NorthEast';
%plt.XTick = -1:0.2:1;
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.ShowBox = 'off';
plt.export(['figures/distance_orderLaw_decay_base1000.pdf']);
%% exponential decay - base 10
% k = 10;
% clear y
% n = 0:1:10;
% for i = 1:length(n)
%     y(i,:) = 100 ./ (1+k.^(d-n(i)));
% end
% figure(1)
% hold on 
% for i = 1:length(n)
%     plot(d,y(i,:),'.-');
% end
% hold off
% xlabel('distance')
% ylabel('dilution level (% of original weight)');
% legend(strsplit(num2str(n)))
% title('decay by base 10')

% %% exponential decay - base 1.5
% k = 1.5;
% clear y
% n = 0:1:10;
% for i = 1:length(n)
%     y(i,:) = 100 ./ (1+k.^(d-n(i)));
% end
% figure(1)
% hold on 
% for i = 1:length(n)
%     plot(d,y(i,:),'.-');
% end
% hold off
% xlabel('distance')
% ylabel('dilution level (% of original weight)');
% legend(strsplit(num2str(n)))
% title('decay by base 1.5')
% %% exponential decay - base 1000
% k = 1000;
% clear y
% n = 0:1:10;
% for i = 1:length(n)
%     y(i,:) = 100 ./ (1+k.^(d-n(i)));
% end
% figure(1)
% hold on 
% for i = 1:length(n)
%     plot(d,y(i,:),'.-');
% end
% hold off
% xlabel('distance')
% ylabel('dilution level (% of original weight)');
% legend(strsplit(num2str(n)))
% title('decay by base 1000')