clear
repeat = 0;
data = [];
while repeat < 1
    
t = 5000;
k = 5; 
A0 =  zeros(5,5,4);
A0(:,:,1) = eye(5);
A0(1,1,1+1) = 0.95*sqrt(2);
A0(1,1,2+1) = -0.9025;
A0(2,1,2+1) = 0.5;
A0(3,1,3+1) = -0.4;
A0(4,1,2+1) = -0.5;
A0(4,4,1+1) = 0.25*sqrt(2);
A0(4,5,1+1) = 0.25*sqrt(2);
A0(5,4,1+1) = -0.25*sqrt(2);
A0(5,5,1+1) = 0.25*sqrt(2);

A = A0;
for i = 2 : 4
    A(:,:,i) = -A0(:,:,i);
end
D = diag([sqrt(0.6),sqrt(0.5),sqrt(0.3),sqrt(0.3),sqrt(0.6)]);
m0 = idarx(A, [], 1);
e = iddata([],  randn(1000+t,k)*D);
temp = sim(m0, e);
data = [data; temp.OutputData(1001:end,:)];
repeat = repeat + 1;
end
data = data';
save('newdata.mat', 'data')
