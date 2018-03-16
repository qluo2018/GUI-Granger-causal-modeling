% input: each column stands for one time series
function postX = preprocessing(X)

% % denoise each column of X
dec = mdwtdec('r',X,2,'db2');
Y = mswden('den',dec,'sqtwolog','sln');
% detrend each column of X
Y = detrend(Y);
% zero mean
Y = cca_rm_temporalmean(Y');
% difference each column of X
X1 = diff(Y');
% X1 = zeros(size(X,1)-1, size(X,2));
% for i = 1 : size(X,2)
%     % denoise each column of X
%     dec = mdwtdec('r',X(:,i),2,'db2');
%     Y = mswden('den',dec,'sqtwolog','sln');
%     % detrend each column of X
%     Y = detrend(Y);
%     % zero mean
%     Y = cca_rm_temporalmean(Y');
%     % difference each column of X
%     Y1 = diff(Y');  
%     X1(:,i) = zscore(Y1);
% end  
% meanX1 = mean(X1,2);
% for i = 1 : size(X,2)
%     X1(:,i) = X1(:,i) - meanX1;
% end
postX =  X1;