function [precision,recall, binaryObs] = prCurve_lump(predictionMat,observationMat,cutoffs,binCallK)
% make the binarized response
binaryObs = [];%up is 1, no change is 0, down is -1
for k= 1:size(observationMat,1)
    observation = observationMat(k,:);
    for i = 1: length(observation)
        for j = i+1: length(observation)
            if observation(j) - observation(i) >= binCallK % up
                binaryObs = [binaryObs,1];
            elseif observation(j) - observation(i) <= -binCallK % down
                binaryObs = [binaryObs,-1];
            else
                binaryObs = [binaryObs,0];
            end
        end
    end
end

% check the prediction
for z = 1:length(cutoffs)
    binaryPred = [];%up is 1, no change is 0, down is -1
    for k = 1: size(predictionMat,1)
        prediction = predictionMat(k,:);
        for i = 1: length(prediction)
            for j = i+1: length(prediction)
                if prediction(j) - prediction(i) >= cutoffs(z) % up
                    binaryPred = [binaryPred,1];
                elseif prediction(j) - prediction(i) <= -cutoffs(z) % down
                    binaryPred = [binaryPred,-1];
                else
                    binaryPred = [binaryPred,0];
                end
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

