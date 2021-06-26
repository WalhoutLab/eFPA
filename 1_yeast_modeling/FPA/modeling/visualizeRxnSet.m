function [Involved_mets1,Dead_ends1] = visualizeRxnSet(rxns,flux,model,exclude)
if strcmp(exclude,'basic')
    excludeMets1 = {'s_0793[ce]','s_0794[c]','s_2783[erm]','s_0795[er]','s_0796[e]','s_3146[gm]','s_0797[g]','s_0798[lp]','s_3094[mm]','s_0799[m]','s_0800[n]','s_0801[p]','s_3164[vm]','s_0802[v]',... H+
        's_3449[ce]','s_0803[c]','s_2808[erm]','s_0804[er]','s_0805[e]','s_2994[gm]','s_0806[g]','s_3657[lp]','s_3226[mm]','s_0807[m]','s_0808[n]','s_0809[p]','s_2976[vm]','s_0810[v]',... h20
        's_1322[c]','s_1323[er]','s_1324[e]','s_1325[g]','s_1326[m]','s_1329[v]','s_3536[ce]','s_2966[erm]','s_2995[gm]','s_2977[vm]','s_3228[mm]',...phosphate
        's_0456[c]','s_2784[erm]','s_0457[er]','s_0458[e]','s_3147[gm]','s_3129[mm]','s_0460[m]','s_0461[n]','s_0462[p]','s_3165[vm]',...co2
    's_2857[ce]','s_0529[c]','s_2785[erm]','s_0530[er]','s_0531[lp]','s_3321[mm]','s_0532[m]','s_0533[n]','s_0534[p]',...coa
    's_0633[c]','s_2834[erm]','s_4157[e]','s_0635[lp]','s_3095[mm]','s_0636[m]','s_0637[n]','s_0638[p]',...ppi
    's_0419[c]','s_0420[e]','s_0421[m]',...nh4
    's_1275[c]','s_2817[erm]','s_1276[er]','s_1277[e]','s_1278[m]','s_1279[p]'...o2
    };
elseif strcmp(exclude, 'all')
    excludeMets1 = {'h[c]','h[m]','h2o[c]','h2o[m]','pi[c]','pi[m]','co2[c]','co2[m]',...
    'coa[c]','coa[m]','ppi[c]','ppi[m]','nh4[c]','o2[c]','o2[m]','atp[c]','atp[m]','adp[c]','adp[m]',...
    'nad[m]','nadp[m]','nadh[m]','nadph[m]','fad[m]','fadh[m]','nad[c]','nadp[c]','nadh[c]','nadph[c]','fad[c]','fadh[c]','sideMet[e]'};
else
    excludeMets1 = {};
end
[Involved_mets1, Dead_ends1] = draw_by_rxn(model, rxns, 'true', 'struc', {''}, excludeMets1, flux);
end
