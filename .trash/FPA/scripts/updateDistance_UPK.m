function [distMat_new,labels_new] = updateDistance_UPK(newUPK,model,distMat,labels,infN)
% update the distance matrix for the newly added UPK reaction. Used in
% the metabolite-centric FPA
%
% the new UPK is assumed to be only allowed for uptaking metabolites
% (regardless of what is set for the boundaries)
%
% USAGE:
%
%    [distMat_new,labels_new] = updateDistance(newDemand,model,distMat,labels,infN)
%
% INPUT:
%    newUPK:           newly added UPK reaction (reaction name) 
%    model:               the COBRA model
%    distMat:             original distance matrix [note: need to be the
%                         raw distance matrix (from rxn A to rxn B,
%                         direction matters)]
%    labels:              labels (the reaction name) of the original
%                         distance matrix
%    infN:                maximum distance (representing infinity)
% OUTPUT:
%    distMat_new:         the new distance matrix 
%    labels_new:          labels for the new distance matrix
%
% .. Author: Xuhang Li, Mar 2020 

% add the new UPK reaction to the end of the distance matrix and label
labels_new = [labels,{[newUPK,'_r']}];
% adding a new dimension to the matrix for the new UPK
distMat_new = zeros(1+size(distMat,1),1+size(distMat,2));
distMat_new(1:size(distMat,1),1:size(distMat,2)) = distMat;
% update the distance of the newUPK [which is min(uptaken met to a
% reaction) + 1]
% first find all reaction that contains the uptaken metabolite
theMet = model.mets(logical(model.S(:,strcmp(model.rxns,newUPK))));
candidateRxns = setdiff(model.rxns(logical(model.S(strcmp(model.mets,theMet),:))),newUPK);
% add direction label
for i = 1:length(candidateRxns)
    if model.S(strcmp(model.mets,theMet),strcmp(model.rxns,candidateRxns{i})) < 0% only the consuming direction is considered
        candidateRxns{i} = [candidateRxns{i},'_f'];
    else
        candidateRxns{i} = [candidateRxns{i},'_r'];
    end
end
% remove the irreversible part (some reactions are actully irreversible)
candidateRxns = intersect(candidateRxns,labels);
% take the minimum and make output
if ~isempty(candidateRxns)
    distMat_new(size(distMat,1)+1,1:size(distMat,2)) = min(distMat(ismember(labels,candidateRxns),:),[],1) + 1;
    distMat_new(1:size(distMat,1),size(distMat,2)+1) = infN; % no reaction can go to upk reaction
else % this UPK is not reachable 
    distMat_new(size(distMat,1)+1,1:size(distMat,2)) = infN;
    distMat_new(1:size(distMat,1),size(distMat,2)+1) = infN; % no reaction can go to upk reaction
end