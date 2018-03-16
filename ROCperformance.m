function [SP, SE, RE, PR] = ROCperformance(score, trueanswer)
for icut = 100 : -1: 1
    TP = 0; TN = 0; FP = 0; FN = 0;
    cutoff = quantile(score,0.01*(icut-1));
    for i = 1 : size(score,1)
        if score(i) > cutoff
            if trueanswer(i) == 1
                TP = TP + 1;
            else
                FP = FP + 1;
            end
        else
            if trueanswer(i) == 0
                TN = TN + 1;
            else
                FN = FN + 1;
            end
        end
    end
    SP(icut) = TN / (TN + FP);
    SE(icut) = TP / (TP + FN);
    RE(icut) = SE(icut);
    PR(icut) = TP / (TP + FP);
end