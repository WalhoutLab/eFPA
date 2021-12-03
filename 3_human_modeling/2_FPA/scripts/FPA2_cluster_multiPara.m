function [FluxPotentials_2d,FluxPotential_solutions_2d] = FPA2_cluster_multiPara(model,targetRxns,master_expression,distMat,labels,nSeq,manualPenalty,manualDist,maxDist,blockList, constantPenalty,parforFlag,penalty_defined,base,alpha)
% Uses the FPA algorithm (`Yilmaz et al., 2020`) to calculate the relative 
% flux potential of a given reaction across conditions. This algorithm 
% finds the objective value of a linear optimization (i.e, maximum flux of 
% a reaction) that best represents the relative expression levels of all
% related gene in certain network neighberhood or the global network. The 
% key concept is to penalizes the flux of reactions according to the 
% relative expression level of those associated genes. A distance order 
% parameter supports to perform such integration at a tunable scale of 
% local metabolic network.
% 
% USAGE:
%
%    FluxPotentials = FPA(model,targetRxns,master_expression,distMat,labels)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    targetRxns:        the target reactions to run FPA on. By default,
%                       both forward and reverse directions of a queried
%                       reaction will be calculated, regardless of the
%                       reaction reversibility. The non-applicable
%                       direction will return NaN in the FPA output.
%    master_expression: the expression profiles of queried conditions. The
%                       master expression variable is required to be a cell
%                       array of structure variables. Each structure
%                       variable corresponds to the expression profile of a
%                       condition in comparison. The structure variable
%                       must have two fields, "genes" and "value"
%                       respectively. For instructions on forming the
%                       master expression variable from a expression table
%                       (i.e, TPM table in text format), please see
%                       "FPA_walkthrough_generic.m"
%    distMat:           the distance matrix for the input model. The matrix
%                       should be distance measures of a irreversible model.
%                       Please see "FPA_walkthrough_generic.m" and 
%                       "MetabolicDistance" section on GitHub on how to 
%                       generate a valid distance matrix
%    labels:            the reaction ID labels for `distMat`. Note that
%                       labels should be for the irreversible version of 
%                       the model. (it will be automatically provided in 
%                       the output of the distance calculator)
%                       
% OPTIONAL INPUTS:
%    n:                 the distance order of FPA calculation. We
%                       suggest 1.5 for C. elegans network. User can vary 
%                       it from 0 (global integration) to 10 (or larger, 
%                       essentially only integrate the expression data 
%                       associated with the target reaction)
%    manualPenalty:     the user-defined penalty for specific reactions.
%                       This input will overide all penalty calculation 
%                       from expression data.
%    manualDist:        the user-defined distance for specific reactions.
%                       This manual distance is define as a single value, 
%                       that is saying, the distance of ALL reactions to 
%                       the specified reaction will be override to the
%                       designated value
%    maxDist:           the user-defined maximum value of the metabolic
%                       distance. All distance values greater than maxDist 
%                       will be overrided with MaxDist. By default, the 
%                       maximum non-infinite value in the distance matrix 
%                       is chosen.
%    blockList:         a list of reactions to block (constrained to zero
%                       flux) during the FPA analysis. Used to conjoin with
%                       iMAT++ result to perform FPA on a context-specific
%                       metabolic network
%    constantPenalty:   A SPECIFIC PARAMETER IN C. ELEGANS DUAL TISSUE
%                       MODELING! It is similar to manualPenalty which 
%                       override the automatically calculated penalties for 
%                       special reactions. However, the constantPenalty 
%                       will NOT override the original penalty of the
%                       super condition, so that FPA of X tissue is 
%                       comparable with that of Intestine. MAY NOT BE
%                       USEFUL GENERAL USE OF FPA.
%     parforFlag:       we support to run the FPA in a parallel manner 
%                       (by default). User can disable the parfor run by
%                       setting the parforFlag to false (we disabled
%                       parfor in metabolite-centric calculation, to
%                       avoid overwhelming time consumption by redundant
%                       penalty calculation). When disable the parfor, user
%                       MUST supply the pre-defined penalty matrix via
%                       "penalty_defined" parameter
%     penalty_defined:  the pre-defined complete penalty matrix for FPA.
%                       Only needed when `parforFlag` is set to false. For
%                       calculating predefined penalty matrix, see
%                       "calculatePenalty.m"
%
%
% OUTPUT:
%   FluxPotentials:     the raw flux potential values of the target
%                       reactions. Potentials are given for both forward 
%                       and reverse direction of each reaction; The column 
%                       order is the same order for the `master_expression`
%                       input (each input conditions). The last column is  
%                       the flux potential of the super condition. For best
%                       evaluation of flux potential, we recommend users to 
%                       normalize the raw flux potential values of each 
%                       condition to the corresponding value of super 
%                       condition, which gives relative flux potential
%                       (rFP) between 0 to 1.
% OPTIONAL OUTPUTS:
%   FluxPotential_solutions:    the FPA solution outputs of each flux
%                               potential objective values. This could be 
%                               used to inspect and understand the flux 
%                               distribution of each flux potential value.
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
    % constant penalty will not be applied to super condition 
    constantPenalty = {};
end
if (nargin < 12)
    % allow non-parallel run for large-scale metabolite-centric optimization;
    % for non-parfor run, the penalty matrix should be pre-defined
    parforFlag = true;
end
if (nargin < 13)
    % allow non-parallel run for large-scale metabolite-centric optimization;
    % for non-parfor run, the penalty calculation should be pre-defined 
    penalty_defined = {};
end
if (nargin < 14)
    % allow non-parallel run for large-scale metabolite-centric optimization;
    % for non-parfor run, the penalty calculation should be pre-defined 
    base = 2;
end
if (nargin < 15)
    alpha = 1;
end

npool = 2;
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

% ONLY FOR C. ELEGANS TISSUE MODELING
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
% unify nomenclature - the nomenclature of reactions should be consistent 
% with the distance calculator
tmp_ind = ~cellfun(@(x) any(regexp(x,'_(f|b|r)$')),model_irrev.rxns); % reactions that have no suffix
model_irrev.rxns(tmp_ind) = cellfun(@(x) [x,'_f'],model_irrev.rxns(tmp_ind), 'UniformOutput',false);% all "_f"
model_irrev.rxns = regexprep(model_irrev.rxns, '_b$','_r');
% create the full distance matrix
distMat2 = maxDist * ones(length(model_irrev.rxns),length(model_irrev.rxns)); % default distance is the maximum distance
[A B] = ismember(model_irrev.rxns,labels);
% filter all Inf values
distMat(distMat > maxDist) = maxDist;
distMat2(A,A) = distMat(B(A),B(A)); % overlay the input distance matrix to default in case any distance is missing in the input
if ~isempty(manualDist)
    [A B] = ismember(model_irrev.rxns,manualDist(:,1));
    distMat2(A,:) = repmat(cell2mat(manualDist(B(A),2)),1,size(distMat2,2)); % overlay the manual-defined distance (fixed distance) onto the matrix
    distMat2(:,A) = repmat(cell2mat(manualDist(B(A),2)),1,size(distMat2,2))'; % overlay the manual-defined distance (fixed distance) onto the matrix
end
% update the distance matrix and label
labels = model_irrev.rxns;
distMat = distMat2;
% reorder expression penalty matrix to match the distance matrix (which
% is ordered to fit in the irreversible model)
penalty_new = ones(length(labels), size(penalty,2));
for i = 1:length(labels)
    myrxn = labels{i};
    myrxn = myrxn(1:(end-2)); % remove the direction suffix
    penalty_new(i,:) = penalty(strcmp(model.rxns,myrxn),:);
end
% update the penalty matrix to the reordered matrix
penalty = penalty_new;
%% part3 merge penalty and calculate the FP
% form and calculate the Flux Potential (which is a linear problem, aka, FBA problem)
fprintf('FPA calculation started:\n');
if parforFlag
    for bInd = 1:length(base)
        k_base = base(bInd);
        for nInd = 1:length(nSeq)
            fprintf(['\n' repmat('.',1,length(targetRxns)) '\n\n']);
            FluxPotentials = cell(length(targetRxns),size(penalty,2));
            FluxPotential_solutions = cell(length(targetRxns),size(penalty,2));
            n = nSeq(nInd);
            environment = getEnvironment();
            parfor i = 1:length(targetRxns)
                restoreEnvironment(environment);
                myrxn = targetRxns{i};
                doForward = any(strcmp(labels, [myrxn,'_f']));% whether to calculate the forward FP, according to the distance matrix  
                doReverse = any(strcmp(labels, [myrxn,'_r']));% whether to calculate the forward FP, according to the distance matrix  

                FPproblems_f = cell(1,size(penalty,2));
                FPproblems_r = cell(1,size(penalty,2));
                for j = 1:size(penalty,2)
                    % block the reactions in the block list
                    model_irrev_tmp0 = model_irrev;
                    if ~isempty(blockList)
                        if j < size(penalty,2)
                            model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
                        else % still block when optimizing supertissue // used to not block, now changed and requires a supercond block list (could be empty)
                            model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
                        end
                    end
                    if doForward 
                        model_irrev_tmp = model_irrev_tmp0;
                        pDist_f = 1./(1+k_base.^(distMat(strcmp(labels, [myrxn,'_f']),:)-n));%the distance term in the weight formula
                        w_f = pDist_f' .* penalty(:,j);% calculate final weight
                        % block the reverse rxns to avoid self-loop
                        model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                        % the FPA is calculated by solvePLP(Penalized Linear Problem) function
                        FPproblems_f{j} = solvePLP_dry(model_irrev_tmp,w_f, labels, [myrxn,'_f'],1, 'max',alpha);
                    else
                        FPproblems_f{j} = {};
                    end
                    if doReverse 
                        model_irrev_tmp = model_irrev_tmp0;
                        pDist_r = 1./(1+k_base.^(distMat(strcmp(labels, [myrxn,'_r']),:)-n));
                        w_r = pDist_r' .* penalty(:,j);
                        % block the forward rxns to avoid self-loop
                        model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                        FPproblems_r{j} = solvePLP_dry(model_irrev_tmp,w_r, labels, [myrxn,'_r'],1, 'max',alpha);
                    else
                        FPproblems_r{j} = {};
                    end
                end

                % solve problems in batch 
                allLP = [FPproblems_f,FPproblems_r];
                allLP_sol = repmat({{}},1,length(allLP));
                allLP_valid = allLP(~cellfun(@isempty,allLP));
                if (~isempty(allLP_valid))
                    batches = 1:npool:length(allLP_valid);
                    if length(batches) == 1
                        solutions = solveGurobi_batch(allLP_valid,'max');
                    else
                        solutions = {};
                        for zz = 1:length(batches)-1
                            solutions = [solutions, solveGurobi_batch(allLP_valid(batches(zz):batches(zz+1)-1),'max')];
                        end
                        solutions = [solutions, solveGurobi_batch(allLP_valid(batches(end):length(allLP_valid)),'max')];
                    end
                    allLP_sol(~cellfun(@isempty,allLP)) = solutions;
                end
                % save result
                FPVector = cell(1,size(penalty,2));
                FPVector_plus = cell(1,size(penalty,2));
                for j = 1:size(penalty,2)
                    if ~isempty(allLP_sol{j})
                        FPVector{1,j}(1) = allLP_sol{j}.obj;
                        % we assume cluster mode is for large dataset, so we dont
                        % save flux distributions 
                        % FPVector_plus{1,j}{1} = allLP_sol{j};
                        % FPVector_plus{1,j}{1}.weight = allLP{j}.weight;
                    else
                        FPVector{1,j}(1) = NaN;
                        % FPVector_plus{1,j}(1) = {''};
                    end
                    if ~isempty(allLP_sol{size(penalty,2) + j})
                        FPVector{1,j}(2) = allLP_sol{size(penalty,2)+j}.obj;
                        % we assume cluster mode is for large dataset, so we dont
                        % save flux distributions 
                        % FPVector_plus{1,j}{2} = allLP_sol{size(penalty,2)+j};
                        % FPVector_plus{1,j}{2}.weight = allLP{size(penalty,2)+j}.weight;
                    else
                        FPVector{1,j}(2) = NaN;
                        % FPVector_plus{1,j}(2) = {''};
                    end
                end

                FluxPotentials(i,:) = FPVector;
                FluxPotential_solutions(i,:) = FPVector_plus;
                fprintf('\b|\n');% for a progress monitor
            end
            if length(nSeq) > 1
                FluxPotentials_all{nInd} = FluxPotentials;
                FluxPotential_solutions_all = {};
            else
                FluxPotentials_all = FluxPotentials;
                FluxPotential_solutions_all = FluxPotential_solutions;
            end
            fprintf('titrating n = %f is done!\n',n);
        end
        if length(base) > 1
            FluxPotentials_2d{bInd} = FluxPotentials_all;
            FluxPotential_solutions_2d = {};
        else
            FluxPotentials_2d = FluxPotentials_all;
            FluxPotential_solutions_2d = FluxPotential_solutions_all;
        end
        fprintf('titrating base = %f is done!\n',k_base);
    end
else % run the same code on a for loop mode
    error('discontinued!')
    for i = 1:length(targetRxns)
        myrxn = targetRxns{i};
        doForward = any(strcmp(labels, [myrxn,'_f']));% whether calculate the forward FP, according to the distance matrix  
        doReverse = any(strcmp(labels, [myrxn,'_r']));
        
        FPproblems_f = cell(1,size(penalty,2));
        FPproblems_r = cell(1,size(penalty,2));
        for j = 1:size(penalty,2)
            % block the reactions in the block list
            model_irrev_tmp0 = model_irrev;
            if ~isempty(blockList)
                if j < size(penalty,2)
                    model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
                else % still block when optimizing supertissue // used to not block, now changed and requires a supercond block list (could be empty)
                    model_irrev_tmp0.ub(ismember(model_irrev_tmp0.rxns,blockList{j})) = 0;
                end
            end
            if doForward 
                model_irrev_tmp = model_irrev_tmp0;
                pDist_f = 1./(1+k_base.^(distMat(strcmp(labels, [myrxn,'_f']),:)-n));
                w_f = pDist_f' .* penalty(:,j);
                % block the reverse rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_r'])) = 0;
                FPproblems_f{j} = solvePLP_dry(model_irrev_tmp,w_f, labels, [myrxn,'_f'],1, 'max',alpha);
            else
                FPproblems_f{j} = {};
            end
            if doReverse 
                model_irrev_tmp = model_irrev_tmp0;
                pDist_r = 1./(1+k_base.^(distMat(strcmp(labels, [myrxn,'_r']),:)-n));
                w_r = pDist_r' .* penalty(:,j);
                % block the forward rxns to avoid self-loop
                model_irrev_tmp.ub(ismember(model_irrev_tmp.rxns,[myrxn,'_f'])) = 0;
                FPproblems_r{j} = solvePLP_dry(model_irrev_tmp,w_r, labels, [myrxn,'_r'],1, 'max',alpha);
            else
                FPproblems_r{j} = {};
            end
        end
        % solve problems in batch 
        allLP = [FPproblems_f,FPproblems_r];
        allLP_sol = repmat({{}},1,length(allLP));
        allLP_valid = allLP(~cellfun(@isempty,allLP));
        if (~isempty(allLP_valid))
            batches = 1:npool:length(allLP_valid);
            if length(batches) == 1
                solutions = solveGurobi_batch(allLP_valid,'max');
            else
                solutions = {};
                for zz = 1:length(batches)-1
                    solutions = [solutions, solveGurobi_batch(allLP_valid(batches(zz):batches(zz+1)-1),'max')];
                end
                solutions = [solutions, solveGurobi_batch(allLP_valid(batches(end):length(allLP_valid)),'max')];
            end
            allLP_sol(~cellfun(@isempty,allLP)) = solutions;
        end
        % save result
        FPVector = cell(1,size(penalty,2));
        FPVector_plus = cell(1,size(penalty,2));
        for j = 1:size(penalty,2)
            if ~isempty(allLP_sol{j})
                FPVector{1,j}(1) = allLP_sol{j}.obj;
                % we assume cluster mode is for large dataset, so we dont
                % save flux distributions 
                % FPVector_plus{1,j}{1} = allLP_sol{j};
                % FPVector_plus{1,j}{1}.weight = allLP{j}.weight;
            else
                FPVector{1,j}(1) = NaN;
                % FPVector_plus{1,j}(1) = {''};
            end
            if ~isempty(allLP_sol{size(penalty,2) + j})
                FPVector{1,j}(2) = allLP_sol{size(penalty,2)+j}.obj;
                % we assume cluster mode is for large dataset, so we dont
                % save flux distributions 
                % FPVector_plus{1,j}{2} = allLP_sol{size(penalty,2)+j};
                % FPVector_plus{1,j}{2}.weight = allLP{size(penalty,2)+j}.weight;
            else
                FPVector{1,j}(2) = NaN;
                % FPVector_plus{1,j}(2) = {''};
            end
        end

        FluxPotentials(i,:) = FPVector;
        FluxPotential_solutions(i,:) = FPVector_plus;
        fprintf('done!\n');
    end
end
