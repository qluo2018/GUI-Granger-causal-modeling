function [timeCau,freq,freqCau,EEout1, EEout2, coeff] = PCauRepeat(timeSeriesX,timeSeriesY,timeSeriesZ,order,sr, Nr, Nl)
%PCAU Summary of this function goes here
%   Detailed explanation goes here
%timeSeriesX Y Z are matrix whose every row is a variable.
%the output is the causality Y to X.
%sr is the sampling rate

[rX,cX] = size(timeSeriesX);
[rY,cY] = size(timeSeriesY);
[rZ,cZ] = size(timeSeriesZ);

[A1,E1] = armorfrepeat([timeSeriesX;timeSeriesZ], Nr, Nl,order);
[A2,E2] = armorfrepeat([timeSeriesX;timeSeriesY;timeSeriesZ],Nr, Nl,order);
EEout = predictionerror(A2, [timeSeriesX;timeSeriesY;timeSeriesZ], Nr, Nl, order);
EEout1 = det(EEout(1:rX, 1:rX));
EEout2 = reshape(E2, 1, (rX+rY+rZ)^2);
coeff = reshape(A2, 1, (rX+rY+rZ)^2*order);
% [E1,A1] = UCM_input([timeSeriesX;timeSeriesZ],order);
% [E2,A2] = UCM_input([timeSeriesX;timeSeriesY;timeSeriesZ],order);

R1 = E1(1:rX,1:rX)-E1(1:rX,rX+1:rX+rZ)*inv(E1(rX+1:rX+rZ,rX+1:rX+rZ))*E1(rX+1:rX+rZ,1:rX);
R2 = E2(1:rX,1:rX)-E2(1:rX,rX+rY+1:rX+rY+rZ)*inv(E2(rX+rY+1:rX+rY+rZ,rX+rY+1:rX+rY+rZ))*E2(rX+rY+1:rX+rY+rZ,1:rX);

timeCau = log(det(R1)/det(R2));

if sr ~= 0
    
%     lamd = 2*pi.*[0:cX/2-1]./cX;  % the original frequency definetion
    lamd = 2*pi.*[0:0.01:sr/2]./(sr/2);     % my definetion of frequency
    
    for m = 1:length(lamd)
        
        %reset matrix for new frequency
        matrixD = eye(rX+rZ);
        matrixB = eye(rX+rY+rZ);
        
        for n = 1:order
            matrixD = matrixD + A1(:,(rX+rZ)*(n-1)+1:(rX+rZ)*n).*exp(-j*n*lamd(m));
            matrixB = matrixB + A2(:,(rX+rY+rZ)*(n-1)+1:(rX+rY+rZ)*n).*exp(-j*n*lamd(m));
        end;
        
        PD = [eye(rX),-E1(1:rX,rX+1:rX+rZ)*inv(E1(rX+1:rX+rZ,rX+1:rX+rZ));zeros(rZ,rX),eye(rZ)];
        %correct matrixD
        matrixD = PD*matrixD;
        
        Exx = E2(1:rX,1:rX);
        Exy = E2(1:rX,rX+1:rX+rY);
        Exz = E2(1:rX,rX+rY+1:rX+rY+rZ);
        Eyy = E2(rX+1:rX+rY,rX+1:rX+rY);
        Eyx = E2(rX+1:rX+rY,1:rX);
        Eyz = E2(rX+1:rX+rY,rX+rY+1:rX+rY+rZ);
        Ezx = E2(rX+rY+1:rX+rY+rZ,1:rX);
        Ezy = E2(rX+rY+1:rX+rY+rZ,rX+1:rX+rY);
        Ezz = E2(rX+rY+1:rX+rY+rZ,rX+rY+1:rX+rY+rZ);
        
        PB1 =   [eye(rX),zeros(rX,rY),-E2(1:rX,rX+rY+1:rX+rY+rZ)*inv(E2(rX+rY+1:rX+rY+rZ,rX+rY+1:rX+rY+rZ));...
            zeros(rY,rX),eye(rY),-E2(rX+1:rX+rY,rX+rY+1:rX+rY+rZ)*inv(E2(rX+rY+1:rX+rY+rZ,rX+rY+1:rX+rY+rZ));...
            zeros(rZ,rX),zeros(rZ,rY),eye(rZ)];
        PB2 =   [eye(rX),zeros(rX,rY),zeros(rX,rZ);...
            -(Exy-Exz*inv(Ezz)*Ezy)*inv(Exx-Exz*inv(Ezz)*Ezx),eye(rY),zeros(rY,rZ);...
            zeros(rZ,rX),zeros(rZ,rY),eye(rZ)];
        PB = PB2*PB1;
        
        matrixB = PB*matrixB;
        matrixH = inv(matrixB);
        
        Qxx = matrixD(1:rX,1:rX)*matrixH(1:rX,1:rX)+matrixD(1:rX,rX+1:rX+rZ)*matrixH(rX+rY+1:rX+rY+rZ,1:rX);
        Qxy = matrixD(1:rX,1:rX)*matrixH(1:rX,rX+1:rX+rY)+matrixD(1:rX,rX+1:rX+rZ)*matrixH(rX+rY+1:rX+rY+rZ,rX+1:rX+rY);
        Qxz = matrixD(1:rX,1:rX)*matrixH(1:rX,rX+rY+1:rX+rY+rZ)+matrixD(1:rX,rX+1:rX+rZ)*matrixH(rX+rY+1:rX+rY+rZ,rX+rY+1:rX+rY+rZ);
        
        %correct error
        cExx = Exx-Exz*inv(Ezz)*Ezx;
        cEyy = Eyy-Eyz*inv(Ezz)*Ezy-(Eyx-Eyz*inv(Ezz)*Ezx)*inv(Exx-Exz*inv(Ezz)*Ezx)*(Exy-Exz*inv(Ezz)*Ezy);
        cEzz = Ezz;
        
        freqCau(m) = log(1+abs(Qxy*cEyy*conj(Qxy')+Qxz*cEzz*conj(Qxz'))/abs(Qxx*cExx*conj(Qxx')));
        
    end;
    
%     freq = sr.*[0:cX/2-1]./cX;     % the original frequency definetion
     freq = 0:0.01:sr/2;
    
end;
