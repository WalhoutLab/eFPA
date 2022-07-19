function model = iCEL2model_xl(fileName, biomassRxnEquation, defaultbound)
%% This function is adapted from original COBRA function "xls2model". Modifications were made to load iCEL1314 model table
%% NOT TESTED FOR ANY OTHER MODEL TABLES
%% Reads a model from Excel spreadsheet.
%
% USAGE:
%
%    model = xls2model(fileName, biomassRxnEquation, defaultbound)
%
% INPUT:
%    fileName:              xls spreadsheet, with one 'KEGG List' and one 'Metabolite List' tab
%
% OPTIONAL INPUTS:
%    biomassRxnEquation:    .xls may have a 255 character limit on each cell,
%                           so pass the biomass reaction separately if it hits this maximum.
%
%    defaultbound:          the deault bound for lower and upper bounds, if
%                           no bounds are specified in the Excel sheet
% OUTPUT:
%    model:                 COBRA Toolbox model
%
% EXAMPLE:
%
%                   'KEGG List' tab headers (case sensitive):
%
%                     * Required:
%
%                       * 'ID':      HEX1
%                       * 'KEGG':          `1 atp[c] + 1 glc-D[c] --> 1 adp[c] + 1 g6p[c] + 1 h[c]`
%
%                     * Optional:
%
%                       * 'GENES':               (3098.3) or (80201.1) or (2645.3) or ...
%                       * 'Name':       Hexokinase
%                       * 'PATHWAY':         Glycolysis
%                       * 'Reversible':        0 (false) or 1 (true)
%                       * 'Lower bound':       0
%                       * 'Upper bound':       1000
%                       * 'Objective':         0/1
%                       * 'Confidence Score':  0,1,2,3,4
%                       * 'EC':         2.7.1.1,2.7.1.2
%                       * 'MachineID':           R000001
%                       * 'Notes':             BIGG also associated with EC 2.7.1.2
%                       * 'References':        PMID:2043117,PMID:7150652,...
%
%                   'Metabolite List' tab: Required headers (case sensitive): (needs to be complete list of metabolites,
%                   i.e., if a metabolite appears in multiple compartments it has to be represented in multiple rows.
%                   IDs need to overlap with use in BIGG List
%
%                     * Required
%
%                       * 'Abbreviation':      glc-D or glc-D[c]
%                     * Optional:
%
%                       * 'Formula' or formula:   C6H12O6
%                       * 'Charge':                       0
%                       * 'Compartment':                  cytosol
%                       * 'Name':                  D-glucose
%                       * 'MachineID':                      C00031
%                       * 'PubChem ID':                   5793
%                       * 'ChEBI ID':                     4167
%                       * 'InChI string':                 InChI=1/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1
%                       * 'SMILES':                       OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O
%                       * 'HMDB ID':                      HMDB00122
%
% NOTE:
%
%    Optional inputs may be required for input on unix machines.
%
% NOTE:
%
%    Find an example Excel sheet at `docs/source/examples/ExcelExample.xlsx`
%
% .. Authors:
%    - Ines Thiele, 01/02/09
%    - Richard Que, 04/27/10, Modified reading of PubChemID and ChEBIID so that if met
%      has multiple IDs, all are passed to model. Confidence Scores
%      PubChemIDs, and ChEBIIDs, are properly passed as cell arrays.
%    - Ronan Fleming, 08/17/10, Support for unix
%    - Hulda S.H., 10/11/10, Modified reading of xls document.
%      Identifies columns by their headers. Added reading of HMDB ID.
%    - Xuhang Li, 03/17/2020. adapted specifically for iCEL1314 table
warning off

% test if Excel is available
excelInstalled = false;
try
    excelObj = actxserver('Excel.Application');
    excelInstalled = true;
    %h.WorkBooks.Item(fileName).Close;
    fprintf(' > Excel is installed.\n\n');
catch ME
    fprintf(' > Excel is not installed; limit of 10000 reactions.\n\n');
end

if ~exist('defaultbound','var')
    defaultbound = 1000;
end

%assumes that one has an xls file with two tabs
if isunix || ~excelInstalled
    [~, Strings, rxnInfo] = xlsread(fileName,'Reaction List', '1:10000');
    [~, MetStrings, metInfo] = xlsread(fileName,'Metabolite List', '1:10000');
    warning on
    if size(Strings, 1) >= 10000 % limit set to 10,000 to prevent out-of-memory issues
        warning('The XLS format is not recommended for large models.')
        warning('Maximum number of reactions reached. Model reaction list truncated at 9,999 reactions.')
    end
    if size(MetStrings, 1) >= 10000 % limit set to 10,000 to prevent out-of-memory issues
        warning('The XLS format is not recommended for large models.')
        warning('Maximum number of metabolites reached. Model metabolite list truncated at 9,999 metabolites.')
    end
    warning off
else
    [~, Strings, rxnInfo] = xlsread(fileName, 'Reaction List');
    [~, MetStrings, metInfo] = xlsread(fileName, 'Metabolite List');
end

%trim empty row from Numbers and MetNumbers
rxnInfo = rxnInfo(1:size(Strings,1),:);
metInfo = metInfo(1:size(MetStrings,1),:);

if isunix && isempty(MetStrings)
    error('Save .xls file as Windows 95 version using gnumeric not openoffice!');
end

requiredRxnHeaders = {'ID','KEGG'};
requiredMetHeaders = {'Abbreviation'};

if ~all(ismember(requiredRxnHeaders,Strings(1,:)))
    error(['Required Headers not present in the "Reaction List" sheet of the provided xls file.', sprintf('\n'),...
           'Note, that headers are case sesnitive!', sprintf('\n'),...
           'Another likely source for this issue is a change in the xls format specification.', sprintf('\n'),...
           'Please have a look at the specification at https://opencobra.github.io/cobratoolbox/docs/ExcelModelFileDefinition.html for the current specifications.']);
end

if ~all(ismember(requiredMetHeaders,MetStrings(1,:)))
    error(['Required Headers not present in the "Metabolite List" sheet of the provided xls file.', sprintf('\n'), ...
           'Note, that headers are case sesnitive!', sprintf('\n'),...
           'Another likely source for this issue is a change in the xls format specification.', sprintf('\n'),...
           'Please have a look at the specification at https://opencobra.github.io/cobratoolbox/docs/ExcelModelFileDefinition.html for the current specifications.']);
end

rxnHeaders = rxnInfo(1,:);

for n = 1:length(rxnHeaders)
    if isnan(rxnHeaders{n})
        rxnHeaders{n} = '';
    end
end

% Assuming first row is header row
rxnAbrList = Strings(2:end,strmatch('ID',rxnHeaders,'exact'));
if ~isempty(strmatch('Name',rxnHeaders,'exact'))
    rxnNameList = Strings(2:end,strmatch('Name',rxnHeaders,'exact'));
else
    rxnNameList = Strings(2:end,strmatch('ID',rxnHeaders,'exact'));
end
rxnList = Strings(2:end,strmatch('KEGG',rxnHeaders,'exact'));
if ~isempty(strmatch('GENES',rxnHeaders,'exact'))
    grRuleList = Strings(2:end,strmatch('GENES',rxnHeaders,'exact'));
else
    grRuleList = cell(size(rxnList,1),1);
    grRuleList(:) = {''};
end

if ~isempty(strmatch('Proteins',rxnHeaders,'exact'))
    Protein = Strings(2:end,strmatch('Proteins',rxnHeaders,'exact'));
end

if ~isempty(strmatch('PATHWAY',rxnHeaders,'exact'))
    subSystemList = Strings(2:end,strmatch('PATHWAY',rxnHeaders,'exact'));
    subSystemList = cellfun(@(x) strsplit(x,';'), subSystemList ,'UniformOutput',0);
else
    subSystemList = cell(size(rxnList,1),1);
    subSystemList(:) = {{''}};
end

if isunix
    for n=1:length(rxnList)
        if length(rxnList{n})==255
            if exist('biomassRxnEquation','var')
                rxnList{n}=biomassRxnEquation;
            else
                error('biomassRxnEquation .xls may have a 255 character limit on each cell, so pass the biomass reaction separately if it hits this maximum.')
            end
        end
    end
end

% initialization with default values
lowerBoundList = -defaultbound*ones(length(rxnAbrList),1);

if ~isempty(strmatch('Lower bound',rxnHeaders,'exact'))
    tmp = rxnInfo(2:end,strmatch('Lower bound',rxnHeaders,'exact'));
    for i = 1:length(tmp)
        if isnumeric(tmp{i})
            lowerBoundList(i) = tmp{i};
        else
            lowerBoundList(i) = str2num(tmp{i});
        end
    end
    lowerBoundList = columnVector(lowerBoundList); %Default -1000
    lowerBoundList(isnan(lowerBoundList)) = -defaultbound;
end

% initialization with default values
upperBoundList = defaultbound*ones(length(rxnAbrList),1);

if ~isempty(strmatch('Upper bound',rxnHeaders,'exact'))
    tmp = rxnInfo(2:end,strmatch('Upper bound',rxnHeaders,'exact'));
    for i = 1:length(tmp)
        if isnumeric(tmp{i})
            upperBoundList(i) = tmp{i};
        else
            upperBoundList(i) = str2num(tmp{i});
        end
    end
    upperBoundList = columnVector(upperBoundList); %Default 1000;
    upperBoundList(isnan(upperBoundList)) = defaultbound;
end

revFlagList = lowerBoundList<0;

% initialization with default values
Objective = zeros(length(rxnAbrList),1);

if ~isempty(strmatch('Objective',rxnHeaders,'exact'))
    tmp = rxnInfo(2:end,strmatch('Objective',rxnHeaders,'exact'));
    for i = 1:length(tmp)
        if isnumeric(tmp{i})
            Objective(i) = tmp{i};
        else
            Objective(i) = str2num(tmp{i});
        end
    end
    Objective = columnVector(Objective);
    Objective(isnan(Objective)) = 0;
end
%boundary information is omited in building initial model. so the
%reversivility is decided by formula
model = createModel_xl(rxnAbrList,rxnNameList,rxnList,'subSystemList',subSystemList,'grRuleList',grRuleList,'lowerBoundList',lowerBoundList,'upperBoundList',upperBoundList);
%%
if ~isempty(strmatch('Confidence Score',rxnHeaders,'exact'))
    model.confidenceScores = cell2mat(rxnInfo(2:end,strmatch('Confidence Score',rxnHeaders,'exact')));
    model.confidenceScores(isnan(model.confidenceScores)) = 0;
end
if ~isempty(strmatch('EC',rxnHeaders,'exact'))    %This needs to be changed to the new annotation scheme and putting the
    %ECNumbers there.
    model.rxnECNumbers = Strings(2:end,strmatch('EC',rxnHeaders,'exact'));
end
if ~isempty(strmatch('Notes',rxnHeaders,'exact'))
    model.rxnNotes = Strings(2:end,strmatch('Notes',rxnHeaders,'exact'));
end
if ~isempty(strmatch('References',rxnHeaders,'exact'))
    model.rxnReferences = Strings(2:end,strmatch('References',rxnHeaders,'exact'));
    numbers = cellfun(@isnumeric ,model.rxnReferences);
    model.rxnReferences(numbers) = cellfun(@convertNumberToID , model.rxnReferences(numbers),'UniformOutput',0);
    model.rxnReferences = cellfun(@(x) regexprep(x,'PMID:',''), model.rxnReferences,'UniformOutput',0);
end

%fill in opt info for metabolites
if ~isempty(Objective) && length(Objective) == length(model.rxns)
    model.c = Objective;
end

metHeaders = metInfo(1,:);

for n = 1:length(metHeaders)
    if isnan(metHeaders{n})
        metHeaders{n} = '';
    end
end

% rename all the metabolite to human readable format
metCol = strmatch('Abbreviation',metHeaders,'exact');
mets = MetStrings(2:end,metCol);
metMachineID = MetStrings(2:end,strmatch('MachineID',metHeaders,'exact'));
metMachineID = cellfun(@(x) x(2:end),metMachineID,'UniformOutput',0);
%lets make a look up table and fill all the annotations
MetInModel = model.mets;
compartSuffix = cellfun(@(x) x(end-2:end),MetInModel,'UniformOutput',0);
MetIDInModel = cellfun(@(x) x(1:min(length(x)-3,regexp(x,'_')-1,'omitnan')),MetInModel,'UniformOutput',0);
TissueLabel = cellfun(@(x) x(regexp(x,'_'):regexp(x,'_')+1),MetInModel,'UniformOutput',0);%if any tissue label
[A,B] = ismember(MetIDInModel,metMachineID);
%change metabolite name (plus compartment)
metsTemp = mets(B(A));
for i =1:length(MetInModel)
    model.mets(i) = {[metsTemp{i},TissueLabel{i},lower(compartSuffix{i})]};
end
%% update other annotations
MetStrings = MetStrings(2:end,:);
%%Set metNames
if ~isempty(strmatch('Name',metHeaders,'exact'))
    model.metNames = columnVector(MetStrings(B(A),strmatch('Name',metHeaders,'exact')));
end
%%Set Formulas
if ~isempty(strmatch('Formula',metHeaders,'exact'))
    model.metFormulas = columnVector(MetStrings(B(A),strmatch('Formula',metHeaders,'exact')));
end
%%Set Charge
if ~isempty(strmatch('Charge',metHeaders,'exact'))
    model.metCharges = cell2mat(metInfo(B(A)+1,strmatch('Charge',metHeaders,'exact')));
end
% set comments
if ~isempty(strmatch('Comments',metHeaders,'exact'))
    model.Comments= columnVector(MetStrings(B(A),strmatch('Comments',metHeaders,'exact')));
end

if ~isempty(strmatch('BiGG',metHeaders,'exact'))
    model.BiGG= columnVector(MetStrings(B(A),strmatch('BiGG',metHeaders,'exact')));
end

if ~isempty(strmatch('MachineID',metHeaders,'exact'))
    model.MetMachineID = columnVector(MetStrings(B(A),strmatch('MachineID',metHeaders,'exact')));
end
if ~isempty(strmatch('Other names',metHeaders,'exact'))
    model.Other_names = columnVector(MetStrings(B(A),strmatch('Other names',metHeaders,'exact')));
end

if ~isempty(strmatch('MW',metHeaders,'exact'))
    model.MolWeight = columnVector(MetStrings(B(A),strmatch('MW',metHeaders,'exact')));
end

if ~isempty(strmatch('Localization',metHeaders,'exact'))
    model.Localization = columnVector(MetStrings(B(A),strmatch('Localization',metHeaders,'exact')));
end

if ~isempty(strmatch('PubChem ID',metHeaders,'exact'))
    %This is a litte trickier, as PubChemIDs are numbers. So we have to
    %load them differently
    model.metPubChemID = columnVector(metInfo(B(A),strmatch('PubChem ID',metHeaders,'exact')));
    numbers = cellfun(@isnumeric ,model.metPubChemID);
    model.metPubChemID(numbers) = cellfun(@convertNumberToID , model.metPubChemID(numbers),'UniformOutput',0);
end
if ~isempty(strmatch('ChEBI ID',metHeaders,'exact'))
    model.metChEBIID  = columnVector(metInfo(B(A),strmatch('ChEBI ID',metHeaders,'exact')));
    numbers = cellfun(@isnumeric ,model.metChEBIID);
    model.metChEBIID(numbers) = cellfun(@convertNumberToID , model.metChEBIID(numbers),'UniformOutput',0);
end

[~,fileName,extension] = fileparts(fileName);

model.description = 'iCEL1314';

warning on
end

function stringNumber = convertNumberToID(number)
if isnan(number)
    stringNumber = '';
else
    stringNumber = num2str(number);
end
end
