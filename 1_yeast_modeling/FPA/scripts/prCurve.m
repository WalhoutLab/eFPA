function [precision,recall,binaryObs] = prCurve(prediction,observation,cutoffs,binCallK)
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
for z = 1:length(cutoffs)
    binaryPred = [];%up is 1, no change is 0, down is -1
    for i = 1: length(prediction)
        for j = i+1: length(prediction)
            if prediction(i) - prediction(j) >= cutoffs(z) % up
                binaryPred = [binaryPred,1];
            elseif prediction(i) - prediction(j) <= -cutoffs(z) % down
                binaryPred = [binaryPred,-1];
            else
                binaryPred = [binaryPred,0];
            end
        end
    end
    precision(z) = (sum((binaryPred == 1 & binaryObs == 1)) ...
                    + sum((binaryPred == -1 & binaryObs == -1))) /...
                    sum((binaryPred == -1 | binaryPred == 1));
    if isnan(precision(z))
        precision(z) = 0;
    end
    recall(z) = (sum((binaryPred == 1 & binaryObs == 1)) ...
            + sum((binaryPred == -1 & binaryObs == -1))) /...
            sum((binaryObs == -1 | binaryObs == 1));     
end
