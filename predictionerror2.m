% prediction error given by conditional AR model 
% Input -------
% data: each row stands for one variable
% A:    coefficient of AR, given by armorfrepeat.m
% order: AR order
% Nr:   # of repeat 
% Nl:   indeces of each trail
% output ------
% the prediction error of X by using Y and Z with coefficient A 
function error = predictionerror2(A, data, Nr, Nl, order)
if order > 1 && size(A,3) == 1
    k = size(data,1);
    temp = A;
    clear A;
    for j = 1 : order
        A(:,:,j) = temp(:,(j-1)*k+1:j*k);       
    end
end
if size(Nl,1) == 1 && size(Nl,2) == 1
    Nl = [1, Nl];
end
error  = [];
for i = 1 : Nr
    for j = Nl(i,1) : Nl(i,1)+order-1
        error = [error, zeros(size(data,1),1)];
    end
    for j = Nl(i,1)+order : Nl(i, 2)        
        temperr = data(:,j);
        for k = 1 : order
            temperr = temperr + A(:,:,k) * data(:,j-k); 
        end
        error = [error, temperr];
    end
end
