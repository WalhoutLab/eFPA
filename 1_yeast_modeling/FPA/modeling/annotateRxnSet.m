function [rxns_annotated, subsystemInfo] = annotateRxnSet(rxns,model, universe)
% annotate the subsystems and calculate the p-value of enrichment 
%rxn        formula     subsystem
%subsystem  p-value ==> now use FDR, 12202020
if (nargin < 3)
    universe = model.rxns;
end
rxns_annotated(:,1) = rxns;
for i = 1:length(rxns)
    rxns_annotated(i,2) = printRxnFormula_XL(model,'rxnAbbrList',rxns(i),'printFlag',0);
    rxns_annotated(i,3) = printGPRForRxns(model,rxns(i),0);
    tmp = model.subSystems{strcmp(model.rxns,rxns{i})};
    rxns_annotated(i,4) = {tmp{1}};
end
rxns_annotated = sortrows(rxns_annotated,4);
% calculate p-value
tmp = model.subSystems(ismember(model.rxns,rxns));
allSubSys = {};
for i = 1:length(tmp)
    for j = 1:length(tmp{i})
        allSubSys = [allSubSys;tmp{i}{j}];
    end
end
subsystemInfo(:,1) = unique(allSubSys);
pvals = [];
for i = 1:length(subsystemInfo(:,1))
    mysys = subsystemInfo(i,1);
    n_obs = 0;
    for j = 1: length(tmp)
        if any(strcmp(tmp{j},mysys{:}))
            n_obs = n_obs +1;
        end
    end
    n_sample = length(rxns);
    n_pop = length(universe);
    n_trueInPop = 0;
    universe_subsys = model.subSystems(ismember(model.rxns,universe));
    for j = 1:length(universe_subsys)
        if any(strcmp(universe_subsys{j},mysys{:}))
            n_trueInPop = n_trueInPop+1;
        end
    end
    subsystemInfo{i,2} = n_obs;
    subsystemInfo{i,3} = n_trueInPop;
    pvals(i) = 1-hygecdf(n_obs-1,n_pop,n_trueInPop,n_sample);
end
fdr = mafdr(pvals,'BHFDR',1);
subsystemInfo(:,4) = mat2cell(fdr',ones(size(subsystemInfo,1),1),1);
subsystemInfo = sortrows(subsystemInfo,4);
end