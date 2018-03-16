function timeS = samplingAR(A,E,origData,Nv,order,Nr,Ns)


Ncache = order;

timeS = [];
for a = 1:Nr
    temptimeS = [];
    for j = 1:order
        temptimeS(:,j) = origData(:,j+(a-1)*Ns);
    end;
    for i = order+1:Ns+Ncache
        pastValue = fliplr(temptimeS(:,i-order:i-1));
        pastValue = pastValue(:);
        temptimeS(:,i) = -A*pastValue+sqrt(diag(E)).*randn(Nv,1);
    end;
    temptimeS(:,1:Ncache) = [];
    timeS = [timeS,temptimeS];
end;
