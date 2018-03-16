function output = wholeTSanalysisSeed(Data, seed)
%%  Whole TS analysis
%%  investigate the connectivity between the seed region and the other
%%  regions
% Input: Data: time*regions*subjects
%        seed: the seed region
% Output: CauREARall (pairwise causality from row to column for each subject)
%         coeffall   (the causal coefficient from row to column)

CauREARall = zeros(size(Data,2),size(Data,2),size(Data,3));
coeffall = zeros(size(Data,2),size(Data,2),size(Data,3));
pvalue = zeros(size(Data,2),size(Data,2),size(Data,3));
order = 1;
for sub = 1 : size(Data,3)
    sub
    % causality row --> column
    for i = 1 : size(seed,2)
        index = seed(i);
        for j = 1 : size(Data,2)
            if j ~= index
                [CauREARall(index,j,sub),r1,r2,coeffall(index,j,sub), pvalue(index,j,sub)]...
                    = Cau(Data(:,j,sub)',Data(:,index,sub)',order);
                [CauREARall(j,index,sub),r1,r2,coeffall(j,index,sub), pvalue(j,index,sub)]...
                    = Cau(Data(:,index,sub)',Data(:,j,sub)',order);
            end
        end
    end
end

output.CauREARall = CauREARall;
output.coeffall = coeffall;
output.pvalue = pvalue;