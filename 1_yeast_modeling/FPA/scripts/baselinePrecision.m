function [precision] = baselinePrecision(observation,binaryCallTol)
% make the binarized response
binaryObs = [];%up is 1, no change is 0, down is -1
for i = 1: length(observation)
    for j = i+1: length(observation)
        if observation(i) - observation(j) >= binaryCallTol(i) % up
            binaryObs = [binaryObs,1];
        elseif observation(i) - observation(j) <= -binaryCallTol(i) % down
            binaryObs = [binaryObs,-1];
        else
            binaryObs = [binaryObs,0];
        end
    end
end
precision = max([sum(binaryObs==1)/length(binaryObs),sum(binaryObs==-1)/length(binaryObs)]);