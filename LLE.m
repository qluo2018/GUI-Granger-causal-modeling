%% fit the AR models with optimal time windows by AIC
function [optcoeff, opterror, optchangepoint, AIC, BIC] = opttvCau(data, Nr, Nl, order, tvCauoptions)
maxpoints = 4;
if nargin <= 4 || isempty(tvCauoptions)
    options.Algorithm = 'interior-point'; %'Active-set';%'trust-region-reflective'; % 
    options.LargeScale = 'on';
    options.Display='off';
    options.Diagnostics='off';
    options.TolX=1e-4;
    options.TolFun=1e-4;
    options.UseParallel = 'Always';
     options.DiffMinChange = 0.002;
else
    options = tvCauoptions;
end

[coeff{1}, error{1}] = armorf(data,Nr,Nl,order);
clear errors
errors = predictionerror2(coeff{1}, data, Nr, [1,Nl], order);
LLF = arma_LLF(errors', error{1}, Nr, [1,Nl], order);
[AIC(1),BIC(1)] = aicbic(LLF,size(data,1)*size(data,1),size(data,2));


for numpoints = 2 : maxpoints
    % initialization
    changepointsini = initializelamenda(numpoints, Nl);    
    % optimization
    parameters = parameterfitting(changepointsini, data, order, Nr, Nl, options);
    % calculate AIC
    num = size(parameters,2);   % number of parameters
    changepoints{numpoints} = recoverchangepoint(parameters, Nl);
    [coeff{numpoints}, error{numpoints}] = tvarmorf(data,Nr,Nl,order,changepoints{numpoints});
    [AIC(numpoints), BIC(numpoints)] = computeAIC(data, Nr, order, numpoints, coeff{numpoints}, error{numpoints}, changepoints{numpoints});
end
[optaic, optNum] = min(BIC);
optchangepoint = changepoints{optNum};
optcoeff = coeff{optNum};
opterror = error{optNum};


function [aicvalue, bicvalue] = computeAIC(X, Nr, order, num, coeff, E1, changepoints)
LLFar = 0;
for i = 1 : num
    clear data
    data = X(:,changepoints(i,1):changepoints(i,2));
    Nl = changepoints(i,2)-changepoints(i,1)+1;
    clear errors
    errors = predictionerror2(coeff{i}, data, Nr, [1,Nl], order);
    LLFar = LLFar + arma_LLF(errors', E1{i}, Nr, [1,Nl], order);
end
r = size(X,1);
[aicvalue, bicvalue] = aicbic(LLFar,num*r*r,size(X,2));


function changepoints = initializelamenda(numpoints, timeserieslength)
windowlength = floor((timeserieslength-30) / numpoints);
if windowlength == 0
    'error in opt_initializelamenda'
end
for i = 1 : numpoints
   changepoints(1,i) = windowlength * i / timeserieslength;
end


function rowpara = reshapeparafortvar(para)
% input: cells of matrices
% output: reshaped these matrices into a row vector
rowpara = [];
for i = 1 : size(para,2)
    rowpara = [rowpara, reshape(para{i},1, size(para{i},1)*size(para{i},2))];
end

function para = recoverparafortvar(rowpara, num, r)
% the reverse function for reshapeparafortvar
for i = 1 : num
    para{i} = reshape(rowpara(1,r*(i-1)+1:r*i), r, r);
end



function parameters = parameterfitting(inipara, data, order, Nr, Nl, options)
ObjectFunction = @(parameters)tvarma_averror(parameters, data, order, Nr, Nl);
num = size(inipara,2);
A = zeros(num-1, num);
for i = 1 : num-1
    A(i,i) = 1;
    A(i,i+1) = -1;    
end

b = -0.05*ones(1, num-1);
lb = zeros(1, num);
lb(1) = 0.05;
ub = ones(1, num);
ub(num) = 1-0.05;
parameters = fmincon(ObjectFunction, inipara, A,b,[],[],lb,ub,[],options);
% parameters = ga(ObjectFunction,num,A,b,[],[],lb,ub);


function averror = tvarma_averror(parameters, X, order, Nr, Nl)
num = size(parameters,2);   % number of parameters
changepoints = recoverchangepoint(parameters, Nl);
[para, error] = tvarmorf(X,Nr,Nl,order,changepoints);
summ = 0;
for i = 1 : num
    summ = summ + det(error{i}) * (changepoints(i,2) - changepoints(i,1)+1);
end
averror = summ / num;
% averror = sum(sum(parameters));


function changepoint = recoverchangepoint(rowpara, len)
num = size(rowpara,2);
rowpara = round(len*rowpara);
for i = 1 : num
    if i == 1 
        changepoint(i,1) = 1;
        changepoint(i,2) = rowpara(1);
    elseif i == num    
        changepoint(i,1) = rowpara(num-1)+1;
        changepoint(i,2) = len;        
    else
        changepoint(i,1) = rowpara(i-1)+1;
        changepoint(i,2) = rowpara(i);                    
    end
end


function [paraini, errorini] = tvarmorf(X,Nr,Nl,order,changepoints)
for i = 1 : size(changepoints,1)
    % fit a model for each time window
    clear data
    data = X(:,changepoints(i,1):changepoints(i,2));
    [A{i}, E{i}] = armorf(data, Nr, changepoints(i,2)-changepoints(i,1)+1, order);
end
paraini = A;
errorini = E;
% paraini = reshapeparafortvar(A);
% errorini = reshapeparafortvar(E);


