%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = updateMetaboliteFormula(model)
%
% Reads the data file and updates the metabolite formula information in the model
%
% William T. Scott, Jr.
% Last Update: 2018-08-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load model
cd ..
model = loadYeastModel;

%Load data:
fid = fopen('../ComplementaryData/modelCuration/Missingmetaboliteformulas.tsv','r');
metaboliteData = textscan(fid,'%s %s %s %s %f32 %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

%Update formulas:
for i = 1:length(metaboliteData{1})
    for j = 1:length(model.mets)
        if strcmp(model.metNames{j},metaboliteData{2}{i})
            model.metFormulas{j} = metaboliteData{3}{i};
        end
    end
end
 
% Save model
saveYeastModel(model)
cd modelCuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
