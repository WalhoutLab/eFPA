function FPA2_met_clusterExecuter_simpleDecay(tmpDir,i,j,batchID)
addpath(genpath('~/cobratoolbox/'));
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
load([tmpDir,'/variables.mat'],'model','master_expression','distMat','distMat_raw','labels','nSeq','maxDist','blockList','constantPenalty','manualDist','penalty','targetRxns','environment','alpha');
restoreEnvironment(environment);
targetExRxns = targetRxns(str2num(i):str2num(j));
exRxns = model.rxns(findExcRxns_XL(model));

for i = 1:length(targetExRxns)
    if contains(targetExRxns{i},'NewMet_') % new demand reaction needed for metabolite in question
        %% first, analyze the production potential
        % it is like a demand reaction except for adding a new reaction
        myMet = targetExRxns{i};
        myMet = myMet(8:end);

        % we skip the blocking as the yeast is cultured in limited defined
        % media, so most metabolite cannot be taken up from media, while
        % essential metabolite cannot be blocked 

%         % find the associated end reactions
%         metInd = cellfun(@(x) ~isempty(regexp(x,['^',regexptranslate('escape',myMet),' \[(\w|\s)*\]$'],'once')),model.metNames);
%         myRxns = model.rxns(any(model.S(metInd,:),1));
%         myRxns = intersect(myRxns,exRxns);

        % block these reactions and optimize 
        model_tmp = model;
        % add the new demand 
        % determine the avaiable compartment --> we use cyotsol as default
%         allCmp = regexp(model.metNames(metInd),'\[(\w|\s)*\]$','match');
%         allCmp = unique([allCmp{:}]);
%         targetCmp = cmpOrder(ismember(cmpOrder,allCmp));
%         targetCmp = targetCmp{1};
        metID = model.mets(strcmp(model.metNames,myMet));
        % add demand rxn
        model_tmp = addReaction(model_tmp,'targetDMN','metaboliteList',metID,'stoichCoeffList',-1,'geneRule', 'NA','lowerBound',0,'printLevel',0);
        targetrxn_fullName = 'targetDMN';
        % the distance of this new demand rxn needs to be calculated; Note
        % that the distance is actually 1+min(distance(any rxn that contains this met)
        [distMat_raw_2,labels_2] = updateDistance(targetrxn_fullName,model_tmp,distMat_raw,labels,maxDist);
        % take the min and block intestine
        distMat_2 = distMat_raw_2;
        for ii = 1:size(distMat_2,1)
            for jj = 1:size(distMat_2,2)
                distMat_2(ii,jj) = min([distMat_raw_2(ii,jj),distMat_raw_2(jj,ii)]);
            end
        end
        % update the penalty mat because of adding reaction
        penalty_defined_2 = [penalty; ones(1,size(penalty,2))];
        [FluxPotential_met_f] = FPA2_cluster_simpleDecay_highIO(model_tmp,{targetrxn_fullName},master_expression,distMat_2,labels_2,nSeq, {},manualDist,maxDist,blockList,constantPenalty,0,penalty_defined_2,alpha);

        %% second, analyze the degradation potential
        model_tmp = model;
        % add the new demand 
        model_tmp = addReaction(model_tmp,'targetUPK','metaboliteList',metID,'stoichCoeffList',-1,'geneRule', 'NA','upperBound',0,'printLevel',0);
        targetrxn_fullName = 'targetUPK';
        % the distance of this new demand rxn needs to be calculated; Note
        % that the distance is actually 1+min(distance(any rxn that contains this met)
        [distMat_raw_2,labels_2] = updateDistance_UPK(targetrxn_fullName,model_tmp,distMat_raw,labels,maxDist);
        % take the min and block intestine
        distMat_2 = distMat_raw_2;
        for ii = 1:size(distMat_2,1)
            for jj = 1:size(distMat_2,2)
                distMat_2(ii,jj) = min([distMat_raw_2(ii,jj),distMat_raw_2(jj,ii)]);
            end
        end
        % update the penalty mat because of adding reaction
        penalty_defined_2 = [penalty; ones(1,size(penalty,2))];
        % block these reactions and optimize 
        % EX reactions in iHuman are all standard EX whose reactant is on
        % the left side. so we block the specific direction
        % block the serection direction to force degradation metabolism of
        % the metabolite of interest
        myMet_tmp = regexprep(myMet,' \[cytoplasm\]','');
        metInd = cellfun(@(x) ~isempty(regexp(x,['^',regexptranslate('escape',myMet_tmp),' \[(\w|\s)*\]$'],'once')),model.metNames);
        myRxns = model.rxns(any(model.S(metInd,:),1));
        myRxns = intersect(myRxns,exRxns);
        model_tmp.ub(ismember(model_tmp.rxns,myRxns)) = 0;
        [FluxPotential_met_r] = FPA2_cluster_simpleDecay_highIO(model_tmp,{targetrxn_fullName},master_expression,distMat_2,labels_2,nSeq, {},manualDist,maxDist,blockList,constantPenalty,0,penalty_defined_2,alpha);
    else % a regular transporter rxns objective
        error('not supported yet')
    end
    %% merge the potentials
    for bInd = 1
        for nInd = 1:length(nSeq)
            FluxPotential_tmp = cell(1,size(penalty,2));
            for z = 1:length(FluxPotential_tmp)
                FluxPotential_tmp{1,z} = [FluxPotential_met_f{bInd}{nInd}{z}(1),FluxPotential_met_r{bInd}{nInd}{z}(2)];
            end 
            FP{bInd}{nInd}(i,:) = FluxPotential_tmp;
        end
    end
end
save([tmpDir,'/FP_',batchID,'.mat'],'FP');

end