function [newModel] = collapseX(model, compID, fluxDistribution)
% This function is designed only for dual-tissue modeling in C.
% elegans. It may have the potential to be modified for similar purpose in other
% models and applications. 
% IT IS NOT DESIRED TO USE IN NON-C. ELEGANS MODEL!
%
% This function collapse a specific compartment in
% the input model according to a flux distribution. The collapsed compartment will be represented as a set of new exchange reactions. The net exchange flux
% between the collapsed compartment and others will be converted to constriants of these exchange reactions, so that the flux burden of the compartment in the input flux distribution will be kept in the 
% resulting new model, yet greatly reduces the model size.
%
% USAGE:
%
%    newModel = collapseX(model, compID, fluxDistribution)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    compID:            The 'id' of the compartment to collapse (i.e., for
%                       cellular compartment, "[c]", the id is "c")
%    fluxDistribution:  the input flux distribution that defines the flux burden
%                       of the collapsed compartment
%
% OUTPUT:
%   newModel:           the collapsed model with flux burden set as
%                       exchange reaction boundaries
%
% `Yilmaz, L. S., Li, X., Nanda, S., Fox, B., Schroeder, F., & Walhout, A. J. (2020). Modeling tissue‐relevant Caenorhabditis elegans metabolism at network, pathway, reaction, and metabolite levels. Molecular Systems Biology, 16(10).
% .. Author: - Xuhang Li, Mar 2020

%% create a new model
newModel = model;
%% collapse metabolite and reactions
% removing the reactions
Xrxns = model.rxns(cellfun(@(x) ~isempty(regexp(x, ['_',compID,'$'], 'once')),model.rxns));
newModel = removeRxns(newModel,Xrxns);
% removing the metabolites
Xmets = model.mets(cellfun(@(x) ~isempty(regexp(x, ['_',compID,'\[.\]$'], 'once')),model.mets));
newModel = removeMetabolites(newModel,Xmets);
%% a special notice on genes
% collapsing genes in X is skipped since reshaping the GR matrix is slow
% one can still do it by uncomment the following codes
% Xgenes = cellfun(@(x) ~isempty(regexp(x, ['_',compID,'$'], 'once')),model.genes);
% newModel =  removeFieldEntriesForType(newModel,Xgenes,'genes',numel(newModel.genes));
%% special notice for dual model
% For dual tissue model, a scar is left in the new model: the side-metabolites constraints (named NonsenseMets in model.mets) for X
% compartement is not deleted, but they are now essentially all zero vector
% in the S matrix. So, it will not influence the modeling.
%% add (new) psudo-exchange rxns
% calculate the flux burden for each metabolite based on transport reactions
% find the X_E metabolites
numericalTol = 0; % tolerance could be larger or equal to solver tolerance; zero gives strigency but may cause solver instability
allExMets = model.mets(cellfun(@(x) ~isempty(regexp(x, '\[e\]$', 'once')),model.mets));
allXmets = model.mets(any(model.S(:,ismember(model.rxns,Xrxns)),2)); % all metablites used in a X reaction
targetMets = intersect(allExMets, allXmets);% target metablite to be proxied by a new exchange reaction
% add reactions
for i = 1:length(targetMets)
    % determine the lb and ub
    % consider the total consuming flux
    associatedRxns = model.rxns(any(model.S(strcmp(model.mets, targetMets{i}),:),1));
    associatedRxns = intersect(associatedRxns,Xrxns);
    demandFlux = -(model.S(strcmp(model.mets, targetMets{i}),ismember(model.rxns,associatedRxns)) * fluxDistribution(ismember(model.rxns,associatedRxns)));
    if demandFlux == 0
        deltaV = 0;
    else
        deltaV = numericalTol;
    end
    % construct the new reaction
    newModel = addReaction(newModel,['EX',model.MetMachineID{strcmp(targetMets{i},model.mets)}(2:end),['_',compID]],...
        'metaboliteList',targetMets(i),...
        'stoichCoeffList',-1,...
        'reversible',1,...
        'subSystem', 'Collapsed compartment X',...
        'lowerBound',demandFlux-deltaV,...
        'upperBound',demandFlux+deltaV,'printLevel',0);
end  

end
