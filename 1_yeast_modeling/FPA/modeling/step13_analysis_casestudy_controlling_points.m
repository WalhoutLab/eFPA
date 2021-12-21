%% About 
% the casestudy of prediction mechanisms of gst synthesis flux
% PS. if we decouple ATP from the two reaction, the predictive power is
% gone; the test is done manually and not shown here
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
GRrate = readtable('./../input/YeastJoshua/originalDataTbl/GRrate.xlsx');
fluxMat_normalized = fluxMat_normalized ./ repmat(GRrate.DR_Actual',size(fluxMat_normalized,1),1);
model = loadYeatModel();
% load expression data 
load('output/normalizedLevels_partialExcluded.mat');
labels2 = regexprep(conditions,'_','-');


%% CLEAR: TCA cycle: control point in the middle
ctrPoint = 'r_1022';
controlledRxns = {'r_0302','r_0280','r_0961','r_0713','r_0505','r_0451','r_0773','r_0658','r_0832','r_0831','r_1022','r_0300'};
plot_func_ctrPoint('TCA',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
%% CLEAR: Pyrimidine metabolism: example of competition (with r_1074) and controling (with others): control point at the begining 
ctrPoint = 'r_0821';
controlledRxns = {'r_1074','r_0307','r_0364','r_1072','r_0820','r_0811','r_0973','r_1045','r_0799','r_0363'};
plot_func_ctrPoint('Pyrimidine_metabolism',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
%% CLEAR: Dihydroorotate production: control point at the end
ctrPoint = 'r_0214';
controlledRxns = {'r_0476','r_1667','r_0214','r_0250','r_0216','r_0958'};
plot_func_ctrPoint('Dihydroorotate_production',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
%% CLEAR: glycolysis: dictation or controling? Besides, glycolysis is coexpressed 
% but r_0959 flux is equal to the glycolysis flux! so may be the flux of
% glycolysis is controled by pyruvate decarboxylase or vice versa 
ctrPoint = 'r_0959';% can also be r_0962 or r_0886 or r_1054 but to a weaker extent 
controlledRxns = {'r_0697','r_1136','r_1054','r_0959','r_2115','r_0892','r_0962','r_0486','r_0893','r_0467','r_0450','r_0886'};
plot_func_ctrPoint('glycolysis',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
%% CLEAR: DEPWTMIKS: coflux of amino acid biosynthesis dictated by 1st step in proline syn
% This is a good example that indicates our analysis cannot tell is it
% controlling or being controlled (correlation)
% these are:
% trp biosyn from chorismate, upsteam of chorismate is not in 
% L-threonine, methionine biosyn from aspartate, but cysteine syn is not in
% isoleucine biosyn from pyruvate
% proline biosyn
% lysine biosyn
% asp-glu conversion
% serine (3-phospho-serine)
% so we may need to reannotate the heatmap
ctrPoint = 'r_0468';
controlledRxns = {'r_0548','r_1041','r_0800','r_0215','r_0219','r_0547','r_0353','r_0663','r_0016','r_0669','r_0468','r_0473',...
    'r_0957','r_0542','r_0566','r_0027','r_0203','r_0202','r_1055','r_1838','r_0988','r_0545','r_0211','r_0727','r_0989','r_0891','r_0891','r_0918'};%r_0887 was excluded because of low absolute corr
plot_func_ctrPoint('AA_DEPWTMIKS',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
% another co-operating modules of AA: LHR
% (1a) upstream branch of leucine syn (from pyruvate) with entire Arginine syn
%     (forms a block) 
% (1b) downstream branch of leucine syn with histidine syn  
% (2) Arginine syn with itself and hist syn, same for hist syn (mostly itself) 

% related: upstream branch of aromatic aa (syn of chorismate) mostly with it
% self, also slightly with Histidine syn
%% UNCLEAR: dNTP: negative correlation! what competition?
ctrPoint = 'r_0886'; % r_0491 or glycolysis; could also be glycolysis, use r_0886 as example
controlledRxns = {'r_0799','r_0973','r_0970','r_1045','r_0364','r_0971'};
plot_func_ctrPoint('dNTP',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
%% UNCLEAR: trehalose_biomass? driven by a few outliers. this is uncertain
% r_0195 or r_1051 (same GPR)
ctrPoint = 'r_1051';
controlledRxns = {'r_0109','r_2141','r_2199','r_2140','r_2197','r_0195','r_1051'};
plot_func_ctrPoint('trehalose_biomass',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
%% UNCLEAR: Sulphate - maybe overfitting
% not correlating with itself, likely overfitting
% the the first step of the L-threonine, methionine biosyn does relate to
% sulphate metabolism (ie homoserine)
ctrPoint = 'r_0215'; 
controlledRxns = {'r_0032','r_0154','r_1026','r_1027','r_0813','r_0883'};
plot_func_ctrPoint('Sulphate_metabolism',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
%% WRONG: Pentose Phosphate : driven by a single point --> may be false positive!
% the downstream reaction r_0889 is not corr but highly co-flux, indicating
% the false positive probability (see the plot
% flux_flux_correlation_ctrPoint_r_0091_targetRxn_r_0889, the outlier
% driving the correlation in other rxns is gone)
ctrPoint = 'r_0091';
controlledRxns = {'r_0990','r_0091','r_0466','r_1050','r_0889'};
plot_func_ctrPoint('PPP',ctrPoint,controlledRxns,normalizedLevel_pro_perPro,valid_rxns_pro_perPro,fluxMat_normalized,rxnLabel,model);
