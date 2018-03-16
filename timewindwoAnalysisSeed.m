function output = timewindwoAnalysisSeed(oData, seed)
%% time window analysis
% Input: Data: time*regions*subjects
%       seed: id's of the seed regions
% Output: avCauREAR, cuCauREAR (i,j,w,sub) (pairwise causality from row to column in
% each subject by each time window set for averaging Granger and cumulative Granger respectively
%         coeff   (the causal coefficient from row to column)

order = 1;
[T,R,S] = size(oData);
avCauREAR = zeros(R, R, S);
cuCauREAR = zeros(R, R, S);
for sub = 1 : S
    sub
    % subject
    
    for iSeed = 1 : size(seed,2)
        index = seed(iSeed);
        for jindex = 1 : R
            if jindex ~= index
                % direction
                clear data
                selected = [index,jindex];
                data = oData(:,selected,sub);
                % optimal window for pairs of time series
                [optcoeff, opterror, optchangepoint, BIC, optexit] = opttvCau(data',1,size(data,1),order, 50);
                disp(optchangepoint)
                alloptCP{sub,iSeed,jindex} = optchangepoint;
                %  allBIC(toymodel) = BIC;
                % calculate Granger causality on optchangepoint for both
                % directions
                combinationAR = [1,2;2 1];
                clear window
                for i = 1 : size(optchangepoint, 1)
                    window{1}{i} = optchangepoint(i,1):optchangepoint(i,2);
                end
                if size(optchangepoint, 1) == 0
                    window{1}{1} = 1 : T;
                end
                % sliding window
                for windowid = 1 : 1
                    clear CauREAR
                    clear R1
                    clear R2
                    for i = 1 : size(combinationAR,1)
                        ratio(1,i) = 0;
                        summ1 = 0;
                        summ2 = 0;
                        for wid = 1 : size(window{windowid},2)
                            clear input_data
                            input_data = data(window{windowid}{wid}, combinationAR(i,:))';    % preparing the input data for algorithm
                            [CauREAR(wid,i), R1{wid,i}, R2{wid,i}, coeff, pvalue] = Cau(input_data(2,:),input_data(1,:),order);
                            summ1 = summ1 + (det(R1{wid,i}));
                            summ2 = summ2 + (det(R2{wid,i}));
                            df1 = order;
                            df2 = window{windowid}{wid}(end)-window{windowid}{wid}(1)+1 - order - order - 1;
                            df22(wid) = df2;
                            ratio(1,i) = ratio(1,i) + (R1{wid,i} - R2{wid,i}) / R2{wid,i} * df2 / df1;
                        end
                        % cumulative GC
                        cuCauREAR(selected(combinationAR(i,1)),selected(combinationAR(i,2)), sub) = log(summ1/summ2);
                        R2window(selected(combinationAR(i,1)),selected(combinationAR(i,2)), sub) =  summ2;
                        df1 = size(window{windowid},1) *order;
                        df2 = T - size(window{windowid},1) * order - size(window{windowid},1) *order - 1;
                        F = (summ1 - summ2) / summ2 * df2 / df1;
                        cupvalue(selected(combinationAR(i,1)),selected(combinationAR(i,2)), sub) = 1-fcdf(F,df1,df2);
                        avCauREAR(selected(combinationAR(i,1)),selected(combinationAR(i,2)), sub) = mean(CauREAR(:,i));% average GC
                    end
                    df1 = order;
                    F = ratio;
                    for i = 1 : 2
                        avpvalue(selected(combinationAR(i,1)),selected(combinationAR(i,2)), sub) = sumF2(F(i), size(window{windowid},1), df1, df22);
                    end
                end
            end
        end
    end
end
output.avCauREAR = avCauREAR;
output.cuCauREAR = cuCauREAR;
output.avpvalue = avpvalue;
output.cupvalue = cupvalue;
output.alloptCP = alloptCP;
