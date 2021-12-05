%% About
% in some cases, the FAA run may fail because the timelimit of the cluster 
% resources is reached. This is a technical accident only specific to our 
% enviuronment, and may not happen every time. We rerun the FAA part for 
% failed ones by 'rerunFailedFAA.sh' when it happens. But when some one is
% reproducing our result, he/she should not encounter any trouble and
% should stop at step3