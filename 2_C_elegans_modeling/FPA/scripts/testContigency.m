function [precision, recall] = testContigency(prediction,observation,cutoff,binCallK)
% make the binarized response
binaryObs = [];%up is 1, no change is 0, down is -1
for i = 1: length(observation)
    for j = i+1: length(observation)
        if observation(i) - observation(j) >= binCallK % up
            binaryObs = [binaryObs,1];
        elseif observation(i) - observation(j) <= -binCallK % down
            binaryObs = [binaryObs,-1];
        else
            binaryObs = [binaryObs,0];
        end
    end
end
% check the prediction
binaryPred = [];%up is 1, no change is 0, down is -1
for i = 1: length(prediction)
    for j = i+1: length(prediction)
        if prediction(i) - prediction(j) >= cutoff % up
            binaryPred = [binaryPred,1];
        elseif prediction(i) - prediction(j) <= -cutoff % down
            binaryPred = [binaryPred,-1];
        else
            binaryPred = [binaryPred,0];
        end
    end
end
% fisher exact test
% tbl = table([sum(binaryObs==1&binaryPred==1);sum(binaryObs==1&binaryPred==-1);sum(binaryObs==1&binaryPred==0)],...
%             [sum(binaryObs==-1&binaryPred==1);sum(binaryObs==-1&binaryPred==-1);sum(binaryObs==-1&binaryPred==0)],...
%             [sum(binaryObs==0&binaryPred==1);sum(binaryObs==0&binaryPred==-1);sum(binaryObs==0&binaryPred==0)],...
%     'VariableNames',{'obs_Up','obs_Down','obs_Unchange'},'RowNames',{'pred_Up','pred_Down','pred_Unchange'});
%p=myfisher33(tbl{:,:});
%[h p] = fishertest(tbl(1:2,1:2),'Tail','right');
precision = (sum(binaryObs==1&binaryPred==1) + sum(binaryObs==-1&binaryPred==-1)) /...
            (sum(binaryPred==1) + sum(binaryPred==-1));
recall = (sum(binaryObs==1&binaryPred==1) + sum(binaryObs==-1&binaryPred==-1)) /...
            (sum(binaryObs==1) + sum(binaryObs==-1));  
if isnan(precision)
    precision = 0;
end
% permutation test
% rng(1126)
% precision_perm = zeros(1000,1);
% recall_perm = zeros(1000,1);
% for t = 1:1000
%     binaryObs_perm = binaryObs(randperm(length(binaryObs)));
%     precision_perm(t) = sum( (binaryObs_perm==1&binaryPred==1) | (binaryObs_perm==-1&binaryPred==-1) ) /...
%             sum( (binaryPred==1) | (binaryPred==-1) );
%     recall_perm(t) = sum( (binaryObs_perm==1&binaryPred==1) | (binaryObs_perm==-1&binaryPred==-1) ) /...
%             sum( (binaryObs_perm==1) | (binaryObs_perm==-1) );   
% end
% precision_perm(isnan(precision_perm)) = 0;
% p_precision = (1+sum(precision_perm >= precision))/1001;
% p_recall = (1+sum(recall_perm >= recall))/1001;
end



