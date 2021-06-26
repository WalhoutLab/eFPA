function type = checkVariability(myexp)
AUCk = 0.3473;
k0 = 0.666;
if min(myexp)>k0
    type = 'Low Variability';
else
    sorted = sort(myexp);
    % scale to 0-1
    sorted = (sorted - min(sorted,[],2)) ./ repmat((max(sorted,[],2) - min(sorted,[],2)),1,25);
    steps = sorted(2:end) - sorted(1:end-1);
    steps_sorted = sort(steps); 
    % reordered curve is
    reordered = [];
    reordered(1) = 0;
    for z = 2:length(sorted)
        reordered(z) = sum(steps_sorted(1:z-1));
    end
    AUC = trapz(1:25,reordered)/12;
    if AUC < AUCk
        type = 'Low Diversity';
    else
        type = 'High Diversity';
    end
end

%%  seperate out the jummping fluxes (for precision)
