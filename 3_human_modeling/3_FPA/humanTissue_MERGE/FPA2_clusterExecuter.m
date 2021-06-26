function FPA2_clusterExecuter(tmpDir,i,j,batchID)
addpath(genpath('~/cobratoolbox/'));
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
load([tmpDir,'/variables.mat'],'model','master_expression','distMat','labels','n','maxDist','blockList','constantPenalty','manualDist','penalty','targetRxns','environment','k_base','alpha');
restoreEnvironment(environment);
[FP,~] = FPA2_cluster(model,targetRxns(str2num(i):str2num(j)),master_expression,distMat,labels,n, {},manualDist,maxDist,blockList,constantPenalty,0,penalty,k_base,alpha);
save([tmpDir,'/FP_',batchID,'.mat'],'FP');
end