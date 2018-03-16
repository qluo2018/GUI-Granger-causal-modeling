function [data, Ht] = arma_bekk_simulate(t,k,kx,ky,parameters,p,q)
% INPUTS:
%     t             - Length of data serie to prouce
%     k             - Dimension of series to produce
%     kx, ky        - Dim for x and y
%     parameters    -- para vector
%     p             - ar_order
%     q             - bekk_order
% 
% OUTPUTS:
%     data          - A t by k matrix of zero mean residuals
%     Ht            - A k x k x t 3 dimension matrix of conditional covariances

t=t+500;
%% Reshape the parameters 
% input: k,kx,ky,p,q,parameters
% output: A, B, const
A = parameters(1 : k*k*p);
C = parameters(k*k*p+1 : k*k*p+kx*(kx+1)/2+ky*(ky+1)/2);
B = parameters(k*k*p+kx*(kx+1)/2+ky*(ky+1)/2+1 : end);
tempA = zeros(k,k,p);
for i=1:p
    tempA(:,:,i) = reshape(A((k*k*(i-1)+1):(k*k*i)),k,k);
end
A = tempA;
Cx = ivech(C(1 : kx*(kx+1)/2));
Cx = tril(Cx);
constx = Cx*Cx';
if ky == 0
    Cy = [];
    consty = [];
else
    Cy = ivech(C(1+kx*(kx+1)/2 : end));
    Cy = tril(Cy);
    consty = Cy*Cy';
end
% const = [constx, zeros(kx,ky); zeros(ky,kx), consty];
tempBx = zeros(kx,k,q);
for i=1:q
    tempBx(:,:,i) = reshape(B((kx*k*(i-1)+1):(kx*k*i)),kx,k);
end
if ky == 0
    BY = [];
else
    tempBy = zeros(ky,k,q);
    for i=1:q
        tempBy(:,:,i) = reshape(B((ky*k*(i-1)+1+kx*k*q):(ky*k*i+kx*k*q)),ky,k);
    end
    BY = tempBy;
end
BX = tempBx;


%% first m data
% use arma model to get the initals 
A0(:,:,1) = eye(k);
for i = 2 : 1+p
    A0(:,:,i) = -A(:,:,i-1)';
end
m0 = idarx(A0, [], 1);
e = iddata([],  randn(1000+t,k));
temp = sim(m0, e);
m = max(p,q);
data = temp.OutputData(500:end,:);

% LHS=eye(k);
% 
% for i=1:p
%     LHS=LHS-A(:,:,i);
% end
% 
% for j=1:q
%     LHS=LHS-B(:,:,j);
% end
% 
% invLHS = LHS^(-1);
% 
% U = invLHS * const * invLHS';
% 
% m = max(p,q);
% data = zeros(m,k);
% Ht = zeros(k,k,m);
% for i = 1 : m
%     data(i,:) = rand(1,k) * U^(0.5);
%     Ht(:,:,i)=U;
% end

%% simulation
for i=m+1:t+m;
    data(i,:) = zeros(1,k);
    for j=1:p
        data(i,:) = data(i,:) + data(i-j,:) * A(:,:,j);
    end    
    Htx(:,:,i) = constx;
    if ky > 0
        Hty(:,:,i) = consty;
    end
    for j=1:q
        Z = (data(i-j,:))'*(data(i-j,:));
        Htx(:,:,i) = Htx(:,:,i) + BX(:,:,j)* Z * BX(:,:,j)';
        if ky > 0 
            Hty(:,:,i) = Hty(:,:,i) + BY(:,:,j)* Z * BY(:,:,j)';
        end        
    end
    if ky > 0
        Ht(:,:,i) = [Htx(:,:,i), zeros(kx,ky); zeros(ky,kx), Hty(:,:,i)];
    else
        Ht(:,:,i) = Htx(:,:,i);
    end
    data(i,:) = data(i,:) + randn(1,k) * (Ht(:,:,i))^(0.5);
end
data=data(m+500:t+m-1,:);
Ht=Ht(:,:,m+500:t+m-1);
%%
% m0 = idarx(A0, [], 1);
% e = iddata([],randn(t,k));
% y = sim(m0, e);
% clear data;
% data = y.OutputData;