function [FluxPotentials,FluxPotential_solutions] = FPA2_oriMERGE_base1000_clusterWrapper(model,targetRxns,master_expression,distMat,labels,n,manualPenalty,manualDist,maxDist,blockList, constantPenalty,parforFlag,penalty_defined,alpha)
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
    n = 1.5;
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
    alpha = 1;
end

k_base = 1000;
jobBatch = 1000;

%% part1 prepare expression matrix and penalty 
% calculate the penalty from expression data
if parforFlag
    fprintf('Mapping the expression levels to penalties...\n');
    penalty = calculatePenalty_singleThread(model,master_expression,manualPenalty);
else
    % in non-parfor mode, we need to provide input penalty (it costs a lot time if
    % repeatedly calculates)
    penalty = penalty_defined;
end

%% start job pooling
fprintf('start job pooling...\n');
rng shuffle
tmpDir = ['tmp_',num2str(fix(rand()*100000))];
mkdir(tmpDir);
environment = getEnvironment();
save([tmpDir,'/variables.mat'],'model','master_expression','distMat','labels','n','manualDist','maxDist','blockList', 'constantPenalty','parforFlag','penalty','targetRxns','environment','k_base','alpha','-v7.3','-nocompression');
% make the target rxn pools
splits = 1:fix(length(targetRxns)/jobBatch)+1:length(targetRxns);
for i = 1:length(splits)-1
    batchID = ['batch_',num2str(i)];
    env_cmd = 'module load gurobi/900 && module load matlab/R2019a && ';
    cmd = ['matlab -nodisplay -nosplash -nojvm -r \"FPA2_clusterExecuter ',tmpDir,' ',num2str(splits(i)),' ',num2str(splits(i+1)-1),' ',batchID,'\"'];
    full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
    bsub_cmd = ['bsub -q short -W 1:00 -n 1 -R rusage[mem=15000] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
    cmd_ready = [bsub_cmd, full_cmd,'"'];
    system(cmd_ready);
    pause(0.25);
end
i = length(splits);
batchID = ['batch_',num2str(i)];
env_cmd = 'module load gurobi/900 && module load matlab/R2019a && ';
cmd = ['matlab -nodisplay -nosplash -nojvm -r \"FPA2_clusterExecuter ',tmpDir,' ',num2str(splits(i)),' ',num2str(length(targetRxns)),' ',batchID,'\"'];
full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
bsub_cmd = ['bsub -q short -W 1:00 -n 1 -R rusage[mem=15000] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
cmd_ready = [bsub_cmd, full_cmd,'"'];
system(cmd_ready);
%% start job monitoring
fprintf('job pooling finished, pausing...\n');
pause(300); % wait for all jobs to be submitted
fprintf('start job monitoring...\n');
runningBatches = {};
for i = 1:length(splits)
    runningBatches{i} = ['batch_',num2str(i)];
end

while ~isempty(runningBatches)
    allFiles = dir([tmpDir,'/*.mat']);
    allFiles = {allFiles.name};
    finished = regexprep(allFiles,'^FP_|.mat.?$','');
    runningBatches = setdiff(runningBatches,finished);
    % find the failed runs
    [~,out] = system('bjobs -q short');
    allTerms = strsplit(out,' ');
    failedBatches = setdiff(runningBatches,allTerms);
    % resubmit the failed batches 
    failedBatches = cellfun(@str2num, regexprep(failedBatches,'batch_',''));
    if any(failedBatches == length(splits)) % if last batch needs to be rerun
        i = length(splits);
        batchID = ['batch_',num2str(i)];
        env_cmd = 'module load gurobi/900 && module load matlab/R2019a && ';
        cmd = ['matlab -nodisplay -nosplash -nojvm -r \"FPA2_clusterExecuter ',tmpDir,' ',num2str(splits(i)),' ',num2str(length(targetRxns)),' ',batchID,'\"'];
        full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
        bsub_cmd = ['bsub -q short -W 1:00 -n 1 -R rusage[mem=15000] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
        cmd_ready = [bsub_cmd, full_cmd,'"'];
        system(cmd_ready);
        fprintf('batch %d was resubmitted...\n',i);
    end
    failedBatches = setdiff(failedBatches,length(splits));
    for j = 1:length(failedBatches)
        i = failedBatches(j);
        batchID = ['batch_',num2str(i)];
        env_cmd = 'module load gurobi/900 && module load matlab/R2019a && ';
        cmd = ['matlab -nodisplay -nosplash -nojvm -r \"FPA2_clusterExecuter ',tmpDir,' ',num2str(splits(i)),' ',num2str(splits(i+1)-1),' ',batchID,'\"'];
        full_cmd = [env_cmd, cmd, ' > ',tmpDir,'/',batchID,'.log && wait'];
        bsub_cmd = ['bsub -q short -W 1:00 -n 1 -R rusage[mem=15000] -e ',tmpDir,'/err_',batchID,'.log ','-J ',batchID,' "'];
        cmd_ready = [bsub_cmd, full_cmd,'"'];
        system(cmd_ready);
        fprintf('batch %d was resubmitted...\n',i);
        pause(0.25);
    end
    fprintf('job monitoring round finished; %d batches remain unfinished; checkpoint interval is 120 sec...\n',length(runningBatches));
    pause(120);
end

%% write out the final result 
FluxPotentials = {};
FluxPotential_solutions = {};
for i = 1:length(splits)
    batchID = ['batch_',num2str(i)];
    load([tmpDir,'/FP_',batchID,'.mat'],'FP');
    FluxPotentials = [FluxPotentials;FP];
end
delete([tmpDir,'/*']);
rmdir(tmpDir);
