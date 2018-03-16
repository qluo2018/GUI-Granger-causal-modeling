%% Detect signal-dependent noise in BOLD signal
% NOTE: this detection is model based, and different models with different orders
% can be selected. Models can be chosen from AR\Polynomial\Fourier 
% by uncommenting the corresponding lines and commenting the lines for other models in subfunction claculateError in
% this file. The default model is AR.
% Input -------
% timeseries: each row stands for one variable
% order: model order, i.e. the maximum time lag
% Nr:   # of repeat 
% Nl:   a matrix with Nr rows 2 columns; the first column specifies the
% starting index of this trial, and the second column for the ending index
% output ------
% the correlation coefficient and pvalue for signal-dependent noise

% reference: 
% Luo Q, Ge T, Grabenhorst F, Feng J, Rolls E. Attention-dependent Modulation of Cortical Taste Circuits Revealed by Granger Causality with Signal-dependent Noise. PLoS Computational Biology, 2013. (accepted)


function [cc0,pp0] = BOLD_SDN_Identify(timeseries, Nr, Nl, order)
temp = Nl;
clear Nl
if Nr == 1
    Nl = [1,temp];
else 
    for i = 1 : Nr
        Nl(i,:) = [temp*(i-1)+1,temp*i];        
    end
end
%% try to find SDN in each region
clear preErr
clear preSig
%%  calculate the error in different ways
[preErr, preSig] = claculateError(timeseries, Nr, Nl, order);

%%  calculate the correlation 
x = [];
y = [];
for i = 1 : size(timeseries,1)
    y_t = [];
    x_t = [];
    for j = 1 : Nr
        y_t = [y_t; preErr(i,Nl(j,1)+order+1:Nl(j,2))];
        x_t = [x_t; timeseries(i,Nl(j,1)+order:Nl(j,2)-1)];
    end
    y = [y, y_t];
    x = [x, x_t];
end
x = mean(x.^2,1);
y = var(y,1);
% compute cc and pvalue
[cc,p111] = corr(x',y');
% fit the line and compute the R-squared
[b,bint,err,rint, stats] = regress(y',[ones(size(x')),x']);

% to better show the data points
newdata = filteroutlier([x;y],5);
x = newdata(1,:);
y = newdata(2,:);
plot(x,y,'.');
hold on
p = polyfit(x,y,1);
stepx = (max(x)-min(x))/100;
x0 = min(x)-stepx : stepx : max(x)+stepx;
f = polyval(p,x0);
plot(x0,f,'r-')
text(quantile(x,0.95),quantile(y,0.95), ['CC=', num2str(cc, '%3.2f'),...
    '\newlinep=', num2str(p111, '%3.2e'), '\newlineR^2=', num2str(stats(1),'%3.2f')])
xlabel('Signal')
ylabel('Noise')
title('X_t = a X_{t-1} + r_t');
% title(['X_t = a X_{t-1} + b X^2_{t-1} + r_t']);
% title('X_t = fourier6(X_{t-1}) + r_t')
cc0 = cc;
pp0 = p111;


% prediction error and prediction given by different models
% NOTE: the code needs to be modified in order to use different models
% Input -------
% timeseries: each row stands for one variable
% order: model order, i.e. the maximum time lag
% Nr:   # of repeat 
% Nl:   a matrix with Nr rows 2 columns; the first column specifies the
% starting index of this trial, and the second column for the ending index
% output ------
% the prediction error (error) and the prediction (pre) established by the selected model 
function [preErr, preSig] = claculateError(timeseries, Nr, Nl, order)
% AR model for each dimension
% n = size(timeseries,1);
% for i = 1 : n
%     clear coeffA1;
%     [coffA1,E1] = armorfrepeat(timeseries(i,:),  Nr, Nl, order);
%     [preErr(i,:),preSig(i,:)] = predictionerror3(coffA1, timeseries(i,:), Nr, Nl, order);
% end

% polynomial
n = size(timeseries,1);
for i = 1 : n
    clear y
    X = [];
    y = timeseries(i,order+1:end)';    
    for j = 1 : order
        X = [X,timeseries(i,j:end-order+j-1)']; %first order
        X = [X,timeseries(i,j:end-order+j-1)'.^2];% second order
    end
    [b,bint,preErr(:,i)] = regress(y,X); 
    preErr(:,i) = y - (X*b);
    preSig(:,i) = X*b;
end
preErr = [zeros(order,n);preErr]';
preSig = [zeros(order,n);preSig]';

% % % fourier
% n = size(timeseries,1);
% for i = 1 : n
%     clear y
%     X = [];
%     y = timeseries(i,order+1:end)';        
%     for j = 1 : order
%         X = [X,timeseries(i,j:end-order+j-1)',timeseries(i,j:end-order+j-1)'.^2]; %first order
%     end
%     tdata = y;
%     preSig(:,i) = zeros(size(y));
%     for j = 1 : size(X,2)
%         f{j} = fit(X(:,j),tdata,'fourier6');
%         tdata = y - f{j}(X(:,j));
%         preSig(:,i) = preSig(:,i) + f{j}(X(:,j));
%     end
%     preErr(:,i) = y - preSig(:,i);
% 
% end
% preErr = [zeros(order,n);preErr]';
% preSig = [zeros(order,n);preSig]';


% prediction error given by AR model 
% Input -------
% data: each row stands for one variable
% A:    coefficient of AR, given by armorfrepeat.m
% order: AR order
% Nr:   # of repeat 
% Nl:   a matrix with Nr rows 2 columns; the first column specifies the
% starting index of this trial, and the second column for the ending index
% output ------
% the prediction error (error) and the prediction (pre) by AR model with coefficient A 
function [error, pre] = predictionerror3(A, data, Nr, Nl, order)
if order > 1 && size(A,3) == 1
    k = size(data,1);
    temp = A;
    clear A;
    for j = 1 : order
        A(:,:,j) = temp(:,(j-1)*k+1:j*k);       
    end
end
error  = [];
pre = [];
for i = 1 : Nr
    for j = Nl(i,1) : Nl(i,1)+order-1
        error = [error, zeros(size(data,1),1)];
        pre = [pre, zeros(size(data,1),1)];
    end
    for j = Nl(i,1)+order : Nl(i, 2)        
        temppre = 0;
        for k = 1 : order
            temppre = temppre - A(:,:,k) * data(:,j-k); 
        end
        pre = [pre, temppre];
        temperr = data(:,j) - temppre;
        error = [error, temperr];
    end
end