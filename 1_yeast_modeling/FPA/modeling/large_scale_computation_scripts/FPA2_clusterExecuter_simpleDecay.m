function FPA2_clusterExecuter_simpleDecay(tmpDir,i,j,batchID)
addpath(genpath('~/cobratoolbox/'));
addpath ./input/
addpath ./../scripts/
addpath ./../scripts/oriMERGE/
addpath ./../../bins/
load([tmpDir,'/variables.mat'],'model','master_expression','distMat','labels','nSeq','maxDist','blockList','constantPenalty','manualDist','penalty','targetRxns','environment','alpha');
restoreEnvironment(environment);
[FP] = FPA2_cluster_simpleDecay_highIO(model,targetRxns(str2num(i):str2num(j)),master_expression,distMat,labels,nSeq, {},manualDist,maxDist,blockList,constantPenalty,0,penalty,alpha);
save([tmpDir,'/FP_',batchID,'.mat'],'FP');
end