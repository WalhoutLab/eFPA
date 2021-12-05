function FPA_clusterExecuter(tmpDir,i,j,batchID)
addpath(genpath('~/cobratoolbox/'));
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
load([tmpDir,'/variables.mat'],'model','master_expression','distMat','labels','n','maxDist','blockList','constantPenalty','penalty','manualDist','targetRxns','alpha','environment');
restoreEnvironment(environment);
[FP,~] = FPA_cluster(model,targetRxns(str2num(i):str2num(j)),master_expression,distMat,labels,n, {},manualDist,maxDist,blockList,constantPenalty,0,penalty,alpha);
save([tmpDir,'/FP_',batchID,'.mat'],'FP');
end