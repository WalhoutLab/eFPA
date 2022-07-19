function preParsedGrRules = preparseGPR_xl(grRules)
% preparse model.grRules before parsing the remaining part
% and transforming model.grRules into model.rules
%
% USAGE:
%
%    preParsedGrRules = preparseGPR(grRules)
%
% INPUT:
%    grRules:           grRules cell or single grRule
%
% OUTPUT:
%    preParsedGrRules:  preparsed grRules cell or single grRule
%
% .. Original Author: -  Laurent Heirendt - December 2017
% .. Modified by: Xuhang Li, Mar 2020 
% .. modified regexp formula to make it suitable for GPR parser

    preParsedGrRules = regexprep(grRules, '[\]\}]',')'); %replace other brackets by parenthesis.
    preParsedGrRules = regexprep(preParsedGrRules, '[\[\{]','('); %replace other brackets by parenthesis.
    preParsedGrRules = regexprep(preParsedGrRules,'([\(])\s*','$1'); %replace all spaces after opening parenthesis
    preParsedGrRules = regexprep(preParsedGrRules,'\s*([\)])','$1'); %replace all spaces before closing paranthesis.
    preParsedGrRules = regexprep(preParsedGrRules, '([\)]\s?|\s)\s*(?i)(and)\s*?(\s?[\(]|\s)\s*', '$1and$3'); %Replace all ands
    preParsedGrRules = regexprep(preParsedGrRules, '([\)]\s?|\s)\s*(?i)(or)\s*?(\s?[\(]|\s)\s*', '$1or$3'); %replace all ors
    preParsedGrRules = regexprep(preParsedGrRules, '[\s]?&[\s]?', ' or '); %introduce spaces around ands
    preParsedGrRules = regexprep(preParsedGrRules, '[\s]?\|[\s]?', ' or '); %introduce spaces around ors.

end