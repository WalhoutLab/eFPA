function [rxns_annotated] = annotateRxnSet2(rxns,model)
% annotate the subsystems and calculate the p-value of enrichment 
%rxn        formula     subsystem
%subsystem  p-value
rxns_annotated(:,1) = rxns;
for i = 1:length(rxns)
    rxns_annotated(i,2) = printRxnFormula_XL(model,'rxnAbbrList',rxns(i),'printFlag',0);
    rxns_annotated(i,3) = model.subSystems(strcmp(model.rxns,rxns{i}));
end
end