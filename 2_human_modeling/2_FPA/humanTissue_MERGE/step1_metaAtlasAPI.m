%% About
% grab metabolite information for mapping the metabolite in HMDB reference
% dataset to the metabolites in the model 
%%
load('input/Human-GEM-RIVEN.mat');
metName = {};
metID = {};
metKEGG = {};
metHMDB = {};
metBiGG = {};
parfor i = 1:length(ihuman.mets)
    [~,msg] = system(['curl -X GET "https://www.metabolicatlas.org/api/Human-GEM/metabolite/',ihuman.mets{i},'/" -H "accept: application/json" -H "X-CSRFToken: N7SxUfN8dixtixZ2a1rFrS1OTQnmPHtwqTpIRwvmFvFHAWPMGXyfA7a7ZtDF9I2K"']);
    metinfo = jsondecode(msg);
    metID{i} = metinfo.id;
    metName{i} = metinfo.name;
    if isfield(metinfo.external_databases,'KEGG')
        metKEGG{i} = metinfo.external_databases.KEGG.id;
    else
        metKEGG{i} = 'NA';
    end
    if isfield(metinfo.external_databases,'HMDB')
        metHMDB{i} = metinfo.external_databases.HMDB.id;
    else
        metHMDB{i} = 'NA';
    end
    if isfield(metinfo.external_databases,'BiGG')
        metBiGG{i} = metinfo.external_databases.BiGG.id;
    else
        metBiGG{i} = 'NA';
    end
    i
end
%%
metTbl = table(metName', metID', metKEGG', metHMDB', metBiGG','VariableNames',{'name','ID','KEGG','HMDB','BiGG'});
writetable(metTbl,'input/iHuman_metTable.csv');
%%
save('metTbl.mat','metTbl');



