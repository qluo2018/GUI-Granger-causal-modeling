%% fit the AR models with optimal time windows by AIC
function [optoptcoeff, optopterror, optoptchangepoint, optoptbic, optoptexit] ...
    = opttvCau(data, Nr, Nl, order, N, maxpoints,tvCauoptions)
if nargin <= 6 || isempty(tvCauoptions)
    options.Algorithm = 'interior-point';%'Active-set';% % 'trust-region-reflective'; %; %;%'trust-region-reflective'; % 
    options.LargeScale = 'off';
    options.Display='iter';
    options.Diagnostics='off';
    options.TolX=1e-4;
    options.TolFun=1e-6;
    options.UseParallel = 'Always';
    options.MaxIter = 1000;
    options.DiffMinChange = 0.002;
    options.stepsize = 0.0001;
else
    options = tvCauoptions;
end

[coeff{1}, error{1}] = armorf(data,Nr,Nl,order);
clear errors
errors = predictionerror2(coeff{1}, data, Nr, [1,Nl], order);
LLF = arma_LLF(errors', error{1}, Nr, [1,Nl], order);
[newAIC(1),newBIC(1)] = aicbic(LLF,size(data,1)*size(data,1)*order,size(data,2));
changepoints{1} = [];
optcoeff{1} = coeff{1};
opterror{1} = error{1};
optexit(1) = 1;

for iTrad = 1 : N
    tradeoff = 0.02*iTrad;
    for numpoints = 2 : maxpoints
        % initialization
        changepointsini = initializelamenda(numpoints, Nl);
        % optimization
        [parameters,EXITFLAG] = parameterfitting(changepointsini, data, order, Nr, Nl, tradeoff, options);
        tempEXITFLAG(numpoints) =  EXITFLAG;
        % calculate AIC
        num = size(parameters,2);   % number of parameters
        changepoints{numpoints} = recoverchangepoint(parameters, Nl);
        output = tvarmorf(data,Nr,Nl,order,changepoints{numpoints});
        coeff{numpoints} = output.A;
        error{numpoints} = output.E;
        aveGC{numpoints} = output.avGC;
        [AIC, BIC] = computeAIC(data, Nr, order, numpoints, ...
            coeff{numpoints}, error{numpoints}, changepoints{numpoints});
        newBIC(numpoints) = BIC;% - tradeoff*log(aveGC{numpoints-1});        
    end    
    [optbic(iTrad), optNum] = min(newBIC);
    opttrade(iTrad) = tradeoff;
    optchangepoint{iTrad} = changepoints{optNum};
    optcoeff{iTrad} = coeff{optNum};
    opterror{iTrad} = error{optNum};
    optexit(iTrad) = tempEXITFLAG(optNum);
end
% plot(opttrade, optbic, '-')
% hold on
[optoptbic, optoptNum] = min(optbic);
% plot(opttrade(optoptNum), optoptbic, 'r*')
% hold off
% ylabel('BIC')
% xlabel('lambda')
optopttrade = opttrade(optoptNum);
optoptchangepoint = optchangepoint{optoptNum};
optoptcoeff = optcoeff{optoptNum};
optopterror = opterror(optoptNum);
optoptexit = optexit(optoptNum);
% pause()



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
[aicvalue, bicvalue] = aicbic(LLFar,num*r*r*order,size(X,2)); 


function changepoints = initializelamenda(numpoints, timeserieslength)
windowlength = floor( (timeserieslength-1) / numpoints);
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



function [parameters,EXITFLAG] = parameterfitting(inipara, data, order, Nr, Nl, tradeoff, options)
ObjectFunction = @(parameters)tvarma_averrorpluscausality(parameters, data, order, Nr, Nl, tradeoff);
num = size(inipara,2);
A = zeros(num-1, num);
for i = 1 : num-1
    A(i,i) = 1;
    A(i,i+1) = -1;    
end
% if Nl < 500
%     scale = 0.1;
% elseif Nl >= 500 &&  Nl < 1000
%     scale = 0.05;
% else
scale = size(data,1) * size(data,1) * order  /Nl;
% end
b = -scale*ones(1, num-1);
lb = zeros(1, num);
lb(1) = scale;
ub = ones(1, num);
ub(num) = 1-scale;
[parameters,FVAL,EXITFLAG] = fmincon(ObjectFunction, inipara, A,b,[],[],lb,ub,[],options);
% options = gaoptimset('InitialPopulation', inipara, 'Display', 'off', 'PopulationSize', 20);
% parameters = ga(ObjectFunction,num,A,b,[],[],lb,ub,[], options);
if size(parameters,2) == 0 
    parameters = inipara;
end
EXITFLAG = 1;

% function averror = tvarma_averror(parameters, X, order, Nr, Nl)
% num = size(parameters,2);   % number of parameters
% changepoints = recoverchangepoint(parameters, Nl);
% [para, error] = tvarmorf(X,Nr,Nl,order,changepoints);
% summ = 0;
% for i = 1 : num
%     summ = summ + det(error{i}) * (changepoints(i,2) - changepoints(i,1)+1);
% end
% averror = summ / num;
% % averror = sum(sum(parameters));

function ErrplusGC = tvarma_averrorpluscausality(parameters, X, order, Nr, Nl, tradeoff)
returnErr = 0;
num = size(parameters,2);   % number of parameters
changepoints = recoverchangepoint(parameters, Nl);
for i = 1 : num
    if changepoints(i,2) - changepoints(i,1)+1 <= 2*size(X,1)
        returnErr = 1;        
    end
end
if returnErr == 0
    output = tvarmorf(X,Nr,Nl,order,changepoints);
    avGC = output.avGC;
    error = output.E;
    summ = 0;
    for i = 1 : num
        summ = summ + det(error{i}) * (changepoints(i,2) - changepoints(i,1)+1);
    end
    averror = summ / num;
    ErrplusGC = averror + tradeoff*1/avGC;
else
    ErrplusGC = 1e50;
end



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


function output = tvarmorf(X,Nr,Nl,order,changepoints)
sumTimeCau = 0;
for i = 1 : size(changepoints,1)
    % fit a model for each time window
    clear data
    data = X(:,changepoints(i,1):changepoints(i,2));
    [A{i}, E{i}] = armorf(data, Nr, changepoints(i,2)-changepoints(i,1)+1, order);
    combination = [1,2; 2,1];    
    for j = 1 : size(combination,1)
        timeSeriesX = data(combination(j,1),:);
        timeSeriesY = data(combination(j,2),:);
        timeCau = Cau(timeSeriesX,timeSeriesY,Nr,changepoints(i,2)-changepoints(i,1)+1,order);
        if timeCau >= 0
            sumTimeCau = sumTimeCau + timeCau;
        end
    end
%     [rss, A{i}] = ls_xx(data', order);
%     E{i} = cov(rss);
end
output.A = A;
output.E = E;
output.avGC = sumTimeCau / size(changepoints,1) / 2;  % average over windows and two directions
% paraini = reshapeparafortvar(A);
% errorini = reshapeparafortvar(E);


