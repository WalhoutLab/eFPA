function MET_FPA_clusterExecuter(tmpDir,i,j,batchID)
addpath(genpath('~/cobratoolbox/'));
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
load([tmpDir,'/variables.mat'],'model','master_expression','distMat','distMat_raw','labels','n','maxDist','blockList','constantPenalty','penalty','manualDist','targetRxns','environment','alpha','helperMet');
restoreEnvironment(environment);
targetExRxns = targetRxns(str2num(i):str2num(j));
FP = cell(length(targetExRxns),size(penalty,2));

% prepare the list of environmental transporter rxns
% the transporter (with env) reactions
metComp = regexp(model.metNames,'\[(\w|\s)*\]$','match');
metComp = [metComp{:}]';
EXmets = strcmp(metComp,'[Extracellular]');
EXinvolvedRxns = model.rxns(any(model.S(EXmets,:)~=0,1));
allCmp_iHumanName = unique(regexprep(model.metNames,' \[(\w|\s)*\]$',''));
allCmp_iHumanName = setdiff(allCmp_iHumanName,{'Side Nutrient in Blood Serum','Major Nutrient in Blood Serum'});
metNames = regexprep(model.metNames,' \[(\w|\s)*\]$','');
TSP = [];
for i = 1:length(allCmp_iHumanName)
    myMet_e = {[allCmp_iHumanName{i},' [Extracellular]']};
    metInd_e = ismember(model.metNames,myMet_e);
    metInd_all = ismember(metNames,allCmp_iHumanName(i));
    metInd_non_e = metInd_all & (~metInd_e);
    myRxns_e = model.rxns(any(model.S(metInd_e,:),1));
    myRxns_non_e = model.rxns(any(model.S(metInd_non_e,:),1));
    % we define the transporter as the reactions that contain the same
    % metabolite in [e] and another compartment (cellular) in the same
    % reaction
    candidate = intersect(myRxns_non_e,myRxns_e);
    % check if is on diff side of the reaction
    if ~isempty(candidate)
        for j = 1:length(candidate) % check if is on diff side of the reaction
            if(sign(model.S(metInd_e,strcmp(model.rxns,candidate(j)))) ~= sign(model.S(metInd_non_e,strcmp(model.rxns,candidate(j)))))
                TSP = union(TSP,candidate(j));
            end
        end
    end
end
% some special transporter will be missed, we add back 
envTspRxns = model.rxns(ismember(model.subSystems,{'Transport reactions'}));
envTspRxns = intersect(envTspRxns,EXinvolvedRxns);
envTspRxns = union(TSP, envTspRxns);

exRxns = model.rxns(findExcRxns_XL(model));
cmpOrder = {'[Cytosol]','[Mitochondria]','[Inner mitochondria]',...
            '[Peroxisome]','[Lysosome]','[Golgi apparatus]',...
            '[Endoplasmic reticulum]','[Nucleus]','[Extracellular]'}; % some met like sucrose and pectin only exist in ex cellular

for i = 1:length(targetExRxns)
    if contains(targetExRxns{i},'NewMet_') % new demand reaction needed for metabolite in question
        %% first, analyze the production potential
        % it is like a demand reaction except for adding a new reaction
        myMet = targetExRxns{i};
        myMet = myMet(8:end);
        % find the associated end reactions
        metInd = cellfun(@(x) ~isempty(regexp(x,['^',regexptranslate('escape',myMet),' \[(\w|\s)*\]$'],'once')),model.metNames);
        myRxns = model.rxns(any(model.S(metInd,:),1));
        myRxns = intersect(myRxns,exRxns);
        % block these reactions and optimize 
        model_tmp = model;
        % add the new demand 
        % determine the avaiable compartment 
        allCmp = regexp(model.metNames(metInd),'\[(\w|\s)*\]$','match');
        allCmp = unique([allCmp{:}]);
        targetCmp = cmpOrder(ismember(cmpOrder,allCmp));
        targetCmp = targetCmp{1};
        metID = model.mets(strcmp(model.metNames,[myMet,' ',targetCmp]));
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
        % EX reactions in iHuman are all standard EX whose reactant is on
        % the left side. so we block the specific direction
        % model_tmp.ub(ismember(model_tmp.rxns,myRxns)) = 0;
        % block uptake 
        model_tmp.lb(ismember(model_tmp.rxns,myRxns)) = 0;
        [FluxPotential_met_f] = FPA_cluster(model_tmp,{targetrxn_fullName},master_expression,distMat_2,labels_2,n, {},manualDist,maxDist,blockList,constantPenalty,false,penalty_defined_2,alpha);   
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
        % block the serection direction
        model_tmp.ub(ismember(model_tmp.rxns,myRxns)) = 0;
        % model_tmp.lb(ismember(model_tmp.rxns,myRxns)) = 0;
        [FluxPotential_met_r] = FPA_cluster(model_tmp,{targetrxn_fullName},master_expression,distMat_2,labels_2,n, {},manualDist,maxDist,blockList,constantPenalty,false,penalty_defined_2,alpha);
        %% merge the potentials
        FluxPotential_tmp = cell(1,size(penalty,2));
        for z = 1:length(FluxPotential_tmp)
            FluxPotential_tmp{1,z} = [FluxPotential_met_f{z}(1),FluxPotential_met_r{z}(2)];
        end       
        FP(i,:) = FluxPotential_tmp;
    else % a regular transporter rxns objective
        % find the metabolite being transported/drained
        myMetSet = setdiff(regexprep(model.metNames(model.S(:,strcmp(model.rxns,[targetExRxns{i}]))<0),' \[(\w|\s)*\]$',''),helperMet);
        %% forward direction potential
        % first obtain the X tissue potential  
        % find the associated end reactions (producing the metabolite)
        % analyze if the forward direction is exporting or importing for EACH
        % mets involved, and block the corresponding rxns
        model_tmp = model;
        targetrxn_fullName = [targetExRxns{i}];
        for k = 1:length(myMetSet) % will be skipped if no myMet was found
            myMet = myMetSet{k};
            % determine if the metabolite is exported or imported
            isExport = ismember({[myMet,' [Extracellular]']},model.metNames(model.S(:,strcmp(model.rxns,targetrxn_fullName))>0)); % if the left side of the reaction is the celluar met
            % find all associated rxns
            metInd = cellfun(@(x) ~isempty(regexp(x,['^',regexptranslate('escape',myMet),' \[(\w|\s)*\]$'],'once')),model.metNames);
            myRxns = model.rxns(any(model.S(metInd,:),1));
            if isExport % metabolite is exported
                % BLOCK UPTAKE DIRECTIONS
                ToBlock = setdiff(intersect(myRxns,envTspRxns),...
                    {targetrxn_fullName});% block import (except for the target TCE itself)
                for y = 1:length(ToBlock)
                    if ismember({[myMet,' [Extracellular]']},model.metNames(model.S(:,strcmp(model.rxns,ToBlock{y}))>0))
                        model_tmp.lb(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    else
                        model_tmp.ub(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    end
                end
            else
                % BLOCK EXPORT DIRECTIONS
                ToBlock = setdiff(intersect(myRxns,envTspRxns),...
                    {targetrxn_fullName});% block import (except for the target TCE itself)
                for y = 1:length(ToBlock)
                    if ismember({[myMet,' [Extracellular]']},model.metNames(model.S(:,strcmp(model.rxns,ToBlock{y}))>0))
                        model_tmp.ub(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    else
                        model_tmp.lb(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    end
                end
            end
        end
        [FluxPotential_met_f] = FPA_cluster(model_tmp,{targetrxn_fullName},master_expression,distMat,labels,n, {},manualDist,maxDist,blockList,constantPenalty,false,penalty,alpha);
        %% reverse direction potential
        model_tmp = model;
        targetrxn_fullName = [targetExRxns{i}];
        for k = 1:length(myMetSet) % will be skipped if no myMet was found
            myMet = myMetSet{k};
            % determine if the metabolite is exported or imported
            isExport = ismember({[myMet,' [Extracellular]']},model.metNames(model.S(:,strcmp(model.rxns,targetrxn_fullName))<0)); % if the right side of the reaction is the celluar met
            % find all associated rxns
            metInd = cellfun(@(x) ~isempty(regexp(x,['^',regexptranslate('escape',myMet),' \[(\w|\s)*\]$'],'once')),model.metNames);
            myRxns = model.rxns(any(model.S(metInd,:),1));
            if isExport % metabolite is exported
                % BLOCK UPTAKE DIRECTIONS
                ToBlock = setdiff(intersect(myRxns,envTspRxns),...
                    {targetrxn_fullName});% block import (except for the target TCE itself)
                for y = 1:length(ToBlock)
                    if ismember({[myMet,' [Extracellular]']},model.metNames(model.S(:,strcmp(model.rxns,ToBlock{y}))>0))
                        model_tmp.lb(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    else
                        model_tmp.ub(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    end
                end
            else
                % BLOCK EXPORT DIRECTIONS
                ToBlock = setdiff(intersect(myRxns,envTspRxns),...
                    {targetrxn_fullName});% block import (except for the target TCE itself)
                for y = 1:length(ToBlock)
                    if ismember({[myMet,' [Extracellular]']},model.metNames(model.S(:,strcmp(model.rxns,ToBlock{y}))>0))
                        model_tmp.ub(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    else
                        model_tmp.lb(strcmp(model_tmp.rxns,ToBlock{y})) = 0;
                    end
                end
            end
        end
        [FluxPotential_met_r] = FPA_cluster(model_tmp,{targetrxn_fullName},master_expression,distMat,labels,n, {},manualDist,maxDist,blockList,constantPenalty,false,penalty,alpha);
        %% merge the potentials
        FluxPotential_tmp = cell(1,size(penalty,2));
        for z = 1:length(FluxPotential_tmp)
            FluxPotential_tmp{1,z} = [FluxPotential_met_f{z}(1),FluxPotential_met_r{z}(2)];
            %FluxPotential_solutions_tmp{1,z} = [FluxPotential_solutions_met_f{z}(1),FluxPotential_solutions_met_r{z}(2)];
        end       
        FP(i,:) = FluxPotential_tmp;
    end
end
save([tmpDir,'/FP_',batchID,'.mat'],'FP');
end