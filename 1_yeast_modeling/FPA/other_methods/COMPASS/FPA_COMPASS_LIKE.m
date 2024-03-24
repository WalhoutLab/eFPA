function [FluxPotentials_all,FluxPotential_solutions_all] = FPA_COMPASS_LIKE(model,targetRxns,master_expression,distMat,labels,nSeq,manualPenalty,manualDist,maxDist,blockList, constantPenalty,parforFlag,penalty_defined)
% Uses the FPA algorithm (`Yilmaz et al., 2020`) to calculate the relative flux
% potential of a given reaction across conditions. This algorithm finds the
% objective value of a linear optimization (i.e, maximum flux of a
% reaction) that best represents the relative expression levels of all
% related gene in a certain network neighberhood or the global network. The key
% concept is to penalizes the flux of reactions according to the relative
% expression level of those associated genes. A distance order parameter
% supports to perform such integration at a tunable scale of metabolic
% neighbers.
% 
% USAGE:
%
%    FluxPotentials = FPA(model,targetRxns,master_expression,distMat,labels,n,manualPenalty,manualDist,maxDist,blockList, constantPenalty,parforFlag,penalty_defined)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    targetRxns:        the target reactions to run FPA on. By default,
%                       both forward and reverse directions of a queried
%                       reaction will be calculated, regardless of the
%                       reaction reversibility. The non-applicable
%                       direction will return NaN in the FPA output.
%    master_expression: the expression profiles of queried conditions. The
%                       master expression variable is reqired to be a cell
%                       array of structure variables. Each structure
%                       variable corresponds to the expression profile of a
%                       condition in comparison. The structure variable
%                       should have two fields, "genes" and "value"
%                       respectively. For instructions on forming the
%                       master expression variable from a expression table
%                       (i.e, TPM table in text format), please see
%                       "FPA_walkthrough_generic.m"
%    distMat:           the distance matrix for the input model. The matrix
%                       should be distance measures of a irreversible model. Please see
%                       "FPA_walkthrough_generic.m" and "MetabolicDistance" section on GitHub
%                       for how to generate a valid distance matrix
%    labels:            the reaction ID labels for distMat. Note that
%                       labels should be for the irreversible version of the model. (it will
%                       be automatically provided in the output of the distance calculator)
% OPTIONAL INPUTS:
%    n:                 the distance order for FPA calculation. We
%                       suggest 1.5 for C. elegans network. User can vary it from 0 (global
%                       integration) to 10 (or larger, essentially only integrate the
%                       expression data associated with the target reaction)
%    manualPenalty:     the user defined penalty for specific reactions.
%                       This input will overide all penalty calculation based on expression
%                       data.
%    manualDist:        the user defined distance for specific reactions.
%                       This manual distance is define as a single value, that is saying, the
%                       distance of ALL reactions to the specified reaction will be override
%                       to the designated value
%    maxDist:           the user-defined maximum value of the metabolic
%                       distance. All distance values greater than maxDist will be overrided
%                       with MaxDist. By default, the maximum non-infinite value in the
%                       distance matrix is chosen.
%    blockList:         a list of reactions to block (constrianed to zero
%                       flux) during the FPA analysis. Used to conjoin with iMAT++ result to
%                       perform FPA on a context-specific metabolic network
%    constantPenalty:   a specifc parameter used in dual tissue integration
%                       in C. elegans tissue modeling. It is similar to manualPenalty which override the automatically calculated penalties for special reactions. 
%                       However, the constantPenalty will NOT override the
%                       original penalty of the super condition, so that FPA of X
%                       tissue is comparable with that of Intestine. May
%                       not be useful for generic using of FPA
%     parforFlag:       we support to run the FPA in a parallel manner (by default). User can disable the parfor run by setting the parforFlag to false 
%                       (we used in metabolite-centric calculation, to
%                       avoid overwhelming time consuming in redundant
%                       penalty calculation). When disable the parfor, user
%                       must supply the pre-defined penalty matrix by
%                       "penalty_defined"
%     penalty_defined:  the pre-defined complete penalty matrix for FPA.
%                       Only needed when parforFlag is set to false. For
%                       calculating predefined penalty matrix, see
%                       "calculatePenalty.m"
%
%
% OUTPUT:
%   FluxPotentials:     the raw flux potential values of the target
%                       reactions. Potentials are given for both forward and reverse direction
%                       of each reaction; The column order is the same order as the
%                       master_expression input, giving the FPA for each input conditions. The
%                       last column is the flux potential of the super condition. For best
%                       evaluation of flux potential, we recommand users to normalize the raw
%                       flux potential values of each condition to the corresponding value of
%                       super condition, so that user will get the relative flux potential
%                       (rFP) between 0 to 1.
% OPTIONAL OUTPUTS:
%   FluxPotential_solutions:    the FPA solution outputs of each flux
%                               potential objective values. This could be used to inspect and
%                               understand the flux distribution of maximum flux potential.
%
% `Yilmaz et al. (2020). Final Tittle and journal.
%
% .. Author: - (COBRA implementation) Xuhang Li, Mar 2020
if (nargin < 6)
    nSeq = 1.5;
end
if (nargin < 7)
    manualPenalty = {};
end
if (nargin < 8)
    manualDist = {};
end
if (nargin < 9)
    maxDist = max(distMat(~isinf(distMat)));
end
if (nargin < 10)
    blockList = {};
end
if (nargin < 11)
    % constant penaly will not be applied to super condition 
    constantPenalty = {};
end
if (nargin < 12)
    % allow non-parallel run for metabolite centric optimization;
    % for non-parfor run, the penalty calculation should be pre-defined to
    % save time
    parforFlag = true;
end
if (nargin < 13)
    % allow non-parallel run for metabolite centric optimization;
    % for non-parfor run, the penalty calculation should be pre-defined to
    % save time
    penalty_defined = {};
end
%% part1 prepare expression matrix and penalty 
% calculate the penalty from expression data
if isempty(penalty_defined)
    fprintf('Mapping the expression levels to penalties...\n');
    penalty = calculatePenalty(model,master_expression,manualPenalty);
else
    % in non-parfor mode, we need to provide input penalty (it costs a lot time if
    % repeatedly calculates)
    penalty = penalty_defined;
end

% apply constant penalty to all conditions except for the super condition
% (at the end of the vector (last column)) 
if ~isempty(constantPenalty)
    [A B] = ismember(model.rxns,constantPenalty(:,1));
    penalty(A,1:end-1) = repmat(cell2mat(constantPenalty(B(A),2)),1,size(penalty,2)-1);
end
%% part2 prepare the distance matrix
% prepare the distance matrix
fprintf('Preparing the distance matrix...\n');
model_irrev = convertToIrreversible(model); % convert to irreversible model
% unify naminclature - the naminclature should be consistant with the
% distance calculator
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); %reactions that have no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);%all "_f"
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
% create the full distance matrix
distMat2 = maxDist * ones(length(model_irrev.rxns),length(model_irrev.rxns)); % default distance is the maximum distance
[A B] = ismember(model_irrev.rxns,labels);
% filter all Inf values
distMat(distMat > maxDist) = maxDist;
distMat2(A,A) = distMat(B(A),B(A)); %overlay the input distance matrix to default in case any distance is missing in the input
if ~isempty(manualDist)
    [A B] = ismember(model_irrev.rxns,manualDist(:,1));
    distMat2(A,:) = repmat(cell2mat(manualDist(B(A),2)),1,size(distMat2,2)); %overlay the manual-defined distance (fixed distance) onto the matrix
    distMat2(:,A) = repmat(cell2mat(manualDist(B(A),2)),1,size(distMat2,2))'; %overlay the manual-defined distance (fixed distance) onto the matrix
end
% update the distance matrix and label
labels = model_irrev.rxns;
distMat = distMat2;
% reorder expression penalty matrix to match the distance matrix (which
% is ordered to fit in the irreversible model)
penalty_new = ones(length(labels), size(penalty,2));
for i = 1:length(labels)
    myrxn = labels{i};
    myrxn = myrxn(1:(end-2)); %remove the direction suffix
    penalty_new(i,:) = penalty(strcmp(model.rxns,myrxn),:);
end
% update the penalty matrix to the reordered matrix
penalty = penalty_new;
%% part3 merge penalty and calculate the FP
%form and calculate the Flux Potential (which is a linear problem, aka, FBA problem)
fprintf('FPA calculation started:\n');
if parforFlag
    for nInd = 1:length(nSeq)
        fprintf(['\n' repmat('.',1,length(targetRxns)) '\n\n']);
        FluxPotentials = cell(length(targetRxns),size(penalty,2));
        FluxPotential_solutions = cell(length(targetRxns),size(penalty,2));
        n = nSeq(nInd);
        environment = getEnvironment();
        parfor i = 1:length(targetRxns)
            restoreEnvironment(environment);
            myrxn = targetRxns{i};
            doForward = any(strcmp(labels, [myrxn,'_f']));%whether to calculate the forward efficiency, according to the distance matrix  
            doReverse = any(strcmp(labels, [myrxn,'_r']));%whether to calculate the forward efficiency, according to the distance matrix  
            efficiencyVector = cell(1,size(penalty,2));
            efficiencyVector_plus = cell(1,size(penalty,2));  
            %% calculate the FP of each condition
            % WHEN WE HAVE DIFFERENT CONSTRAINTS IN DIFFERENT CONDITIONS,
            % WE NEED TO CALCULATE THE MINALLOWANCE FOR EACH CONDITION
            % INDIVIDUALLY;  
            for j = 1:(size(penalty,2)) %loop over all the queried conditions
                % block the reactions in the block list
                model_irrev_tmp0 = model_irrev;
                if ~isempty(blockList)
                    if j < size(penalty,2)
                        model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
                    else % still block when optimizing supertissue // used to not block, now changed and requires a supercond block list (could be empty)
                        model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
                    end
                end
                
                % calculate the maximum flux for the target reaction
                if doForward %calculate the minAllowance of the forward reaction              
                    model_irrev_tmp = model_irrev_tmp0;
                    % block the reverse rxns to avoid self-loop
                    model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                    model_irrev_tmp = changeObjective(model_irrev_tmp,[myrxn,'_f']);
                    solVm_f = optimizeCbModel(model_irrev_tmp,'max');
                else
                    solVm_f = NaN;
                end
                if doReverse 
                    model_irrev_tmp = model_irrev_tmp0;
                    % block the forward rxns to avoid self-loop
                    model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                    model_irrev_tmp = changeObjective(model_irrev_tmp,[myrxn,'_r']);
                    solVm_r = optimizeCbModel(model_irrev_tmp,'max');
                else
                    solVm_r = NaN;
                end     
                
                % calculate the real FP
                if doForward && isstruct(solVm_f)                
                    model_irrev_tmp = model_irrev_tmp0;
                    pDist_f = 1./((1+distMat(strcmp(labels, [myrxn,'_f']),:)).^n);%the distance term in the weight formula
                    w_f = pDist_f' .* penalty(:,j);%calculate final weight
                    % block the reverse rxns to avoid self-loop
                    model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                    % the FPA is calculated by solvePLP(Penalized Linear Problem) function
                    [efficiencyVector{1,j}(1),efficiencyVector_plus{1,j}{1}] = GetMinAllowance_COMPASS(model_irrev_tmp,w_f, labels, [myrxn,'_f'],1,solVm_f);
                else
                    efficiencyVector{1,j}(1) = NaN;
                    efficiencyVector_plus{1,j}(1) = {''};
                end
                if doReverse && isstruct(solVm_r)
                    model_irrev_tmp = model_irrev_tmp0;
                    pDist_r = 1./((1+distMat(strcmp(labels, [myrxn,'_r']),:)).^n);
                    w_r = pDist_r' .* penalty(:,j);
                    % block the forward rxns to avoid self-loop
                    model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                    [efficiencyVector{1,j}(2),efficiencyVector_plus{1,j}{2}] = GetMinAllowance_COMPASS(model_irrev_tmp,w_r, labels, [myrxn,'_r'],1,solVm_r);
                else
                    efficiencyVector{1,j}(2) = NaN;
                    efficiencyVector_plus{1,j}(2) = {''};
                end
            end
            FluxPotentials(i,:) = efficiencyVector;
            FluxPotential_solutions(i,:) = efficiencyVector_plus;
            fprintf('\b|\n');%for simple progress monitor
        end
        if length(nSeq) > 1
            FluxPotentials_all{nInd} = FluxPotentials;
            FluxPotential_solutions_all{nInd} = {};
        else
            FluxPotentials_all = FluxPotentials;
            FluxPotential_solutions_all = FluxPotential_solutions;
        end
        fprintf('titrating n = %f is done!\n',n);
    end
else %run the same code on a for loop mode
    
    % discontinued!!!!!!!!!
    
    
    for nInd = 1:length(nSeq)
        FluxPotentials = cell(length(targetRxns),size(penalty,2));
        FluxPotential_solutions = cell(length(targetRxns),size(penalty,2));
        n = nSeq(nInd);
        for i = 1:length(targetRxns)
            myrxn = targetRxns{i};
            doForward = any(strcmp(labels, [myrxn,'_f']));%whether to calculate the forward efficiency, according to the distance matrix  
            doReverse = any(strcmp(labels, [myrxn,'_r']));%whether to calculate the forward efficiency, according to the distance matrix  
            efficiencyVector = cell(1,size(penalty,2));
            efficiencyVector_plus = cell(1,size(penalty,2));
            %% calculate the minimal flux allowance for the target reaction (super condition)
            % block the reactions in the block list
            model_irrev_tmp0 = model_irrev;
            if ~isempty(blockList)
                % apply the block list for the super condition(could be empty)
                model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{size(penalty,2)})) = 0;
            end
            if doForward %calculate the minAllowance of the forward reaction              
                model_irrev_tmp = model_irrev_tmp0;
                pDist_f = 1./((1+distMat(strcmp(labels, [myrxn,'_f']),:)).^n);%the distance term in the weight formula
                w_f = pDist_f' .* penalty(:,size(penalty,2));%calculate final weight
                % block the reverse rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                [minAllowance_f,efficiencyVector{1,size(penalty,2)}(1),efficiencyVector_plus{1,size(penalty,2)}{1}] = GetMinAllowance(model_irrev_tmp,w_f, labels, [myrxn,'_f'],1, 'max');
                % the solution for minimal allowance is also the solution for
                % super condition
            else
                minAllowance_f = NaN;
                efficiencyVector{1,size(penalty,2)}(1) = NaN;
                efficiencyVector_plus{1,size(penalty,2)}(1) = {''};
            end
            if doReverse 
                model_irrev_tmp = model_irrev_tmp0;
                pDist_r = 1./((1+distMat(strcmp(labels, [myrxn,'_r']),:)).^n);
                w_r = pDist_r' .* penalty(:,size(penalty,2));
                % block the forward rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                [minAllowance_r,efficiencyVector{1,size(penalty,2)}(2),efficiencyVector_plus{1,size(penalty,2)}{2}] = GetMinAllowance(model_irrev_tmp,w_r, labels, [myrxn,'_r'],1, 'max');
            else
                minAllowance_r = NaN;
                efficiencyVector{1,size(penalty,2)}(2) = NaN;
                efficiencyVector_plus{1,size(penalty,2)}(2) = {''};
            end       
            %% calculate the FP of each condition
            for j = 1:(size(penalty,2)-1) %loop over all the queried conditions (except for super cond, which we calculated above)
                % block the reactions in the block list
                model_irrev_tmp0 = model_irrev;
                if ~isempty(blockList)
                    model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
                end
                if doForward && ~isnan(minAllowance_f)                
                    model_irrev_tmp = model_irrev_tmp0;
                    pDist_f = 1./((1+distMat(strcmp(labels, [myrxn,'_f']),:)).^n);%the distance term in the weight formula
                    w_f = pDist_f' .* penalty(:,j);%calculate final weight
                    % block the reverse rxns to avoid self-loop
                    model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                    % the FPA is calculated by solvePLP(Penalized Linear Problem) function
                    [efficiencyVector{1,j}(1),efficiencyVector_plus{1,j}{1}] = solvePLP_minAllowance(model_irrev_tmp,w_f, labels, [myrxn,'_f'],1,minAllowance_f, 'max');
                else
                    efficiencyVector{1,j}(1) = NaN;
                    efficiencyVector_plus{1,j}(1) = {''};
                end
                if doReverse && ~isnan(minAllowance_r)
                    model_irrev_tmp = model_irrev_tmp0;
                    pDist_r = 1./((1+distMat(strcmp(labels, [myrxn,'_r']),:)).^n);
                    w_r = pDist_r' .* penalty(:,j);
                    % block the forward rxns to avoid self-loop
                    model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                    [efficiencyVector{1,j}(2),efficiencyVector_plus{1,j}{2}] = solvePLP_minAllowance(model_irrev_tmp,w_r, labels, [myrxn,'_r'],1,minAllowance_r, 'max');
                else
                    efficiencyVector{1,j}(2) = NaN;
                    efficiencyVector_plus{1,j}(2) = {''};
                end
            end
            FluxPotentials(i,:) = efficiencyVector;
            FluxPotential_solutions(i,:) = efficiencyVector_plus;
        end
        if length(nSeq) > 1
            FluxPotentials_all{nInd} = FluxPotentials;
            FluxPotential_solutions_all{nInd} = FluxPotential_solutions;
        else
            FluxPotentials_all = FluxPotentials;
            FluxPotential_solutions_all = FluxPotential_solutions;
        end
        fprintf('titrating n = %f is done!\n',n);
    end
    fprintf('done!\n');%for simple progress monitor
end
