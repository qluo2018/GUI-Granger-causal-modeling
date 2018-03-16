function [timeCau,freq,freqCau] = CCau(timeSeriesX,timeSeriesY,timeSeriesZ,Nr,Nl,order,sr)
%CCAU Summary of this function goes here
%   Detailed explanation goes here
%timeSeriesX Y Z are matrix whose every row is a variable.
%the output is the causality Y to X.
%sr is the sampling rate

[rX,cX] = size(timeSeriesX);
[rY,cY] = size(timeSeriesY);
[rZ,cZ] = size(timeSeriesZ);

[A1,E1] = armorf([timeSeriesX;timeSeriesZ],Nr,Nl,order);
[A2,E2] = armorf([timeSeriesX;timeSeriesY;timeSeriesZ],Nr,Nl,order);
% [E1,A1] = UCM_input([timeSeriesX;timeSeriesZ],order);
% [E2,A2] = UCM_input([timeSeriesX;timeSeriesY;timeSeriesZ],order);

timeCau = log(trace(E1(1:rX,1:rX))/trace(E2(1:rX,1:rX)));

if sr ~= 0
    lamd = 2*pi.*[0:cX/2-1]./cX;

    for m = 1:length(lamd)

        %reset matrix for new frequency
        matrixD = eye(rX+rZ);
        matrixB = eye(rX+rY+rZ);

        for n = 1:order
            matrixD = matrixD + A1(:,(rX+rZ)*(n-1)+1:(rX+rZ)*n).*exp(-j*n*lamd(m));
            matrixB = matrixB + A2(:,(rX+rY+rZ)*(n-1)+1:(rX+rY+rZ)*n).*exp(-j*n*lamd(m));
        end;

        PD = [eye(rX),zeros(rX,rZ);-E1(rX+1:rX+rZ,1:rX)*inv(E1(1:rX,1:rX)),eye(rZ)];
        %correct matrixD
        matrixD = PD*matrixD;

        PB1 =   [eye(rX),zeros(rX,rY),zeros(rX,rZ);...
                -E2(rX+1:rX+rY,1:rX)*inv(E2(1:rX,1:rX)),eye(rY),zeros(rY,rZ);...
                -E2(rX+rY+1:rX+rY+rZ,1:rX)*inv(E2(1:rX,1:rX)),zeros(rZ,rY),eye(rZ)];
        PB2 =   [eye(rX),zeros(rX,rY),zeros(rX,rZ);...
                zeros(rY,rX),eye(rY),zeros(rY,rZ);...
                zeros(rZ,rX),-(E2(rX+rY+1:rX+rY+rZ,rX+1:rX+rY)-E2(rX+rY+1:rX+rY+rZ,1:rX)*inv(E2(1:rX,1:rX))*E2(1:rX,rX+1:rX+rY))*...
                inv(E2(rX+1:rX+rY,rX+1:rX+rY)-E2(rX+1:rX+rY,1:rX)*inv(E2(1:rX,1:rX))*E2(1:rX,rX+1:rX+rY)),eye(rZ)];
        PB = PB2*PB1;

        matrixB = PB*matrixB;

        matrixH = inv(matrixB);

        Qxx = matrixD(1:rX,1:rX)*matrixH(1:rX,1:rX)+matrixD(1:rX,rX+1:rX+rZ)*matrixH(rX+rY+1:rX+rY+rZ,1:rX);
        Qxy = matrixD(1:rX,1:rX)*matrixH(1:rX,rX+1:rX+rY)+matrixD(1:rX,rX+1:rX+rZ)*matrixH(rX+rY+1:rX+rY+rZ,rX+1:rX+rY);
        Qxz = matrixD(1:rX,1:rX)*matrixH(1:rX,rX+rY+1:rX+rY+rZ)+matrixD(1:rX,rX+1:rX+rZ)*matrixH(rX+rY+1:rX+rY+rZ,rX+rY+1:rX+rY+rZ);

        Exx = E2(1:rX,1:rX);
        Exy = E2(1:rX,rX+1:rX+rY);
        Exz = E2(1:rX,rX+rY+1:rX+rY+rZ);
        Eyy = E2(rX+1:rX+rY,rX+1:rX+rY);
        Eyx = E2(rX+1:rX+rY,1:rX);
        Eyz = E2(rX+1:rX+rY,rX+rY+1:rX+rY+rZ);
        Ezx = E2(rX+rY+1:rX+rY+rZ,1:rX);
        Ezy = E2(rX+rY+1:rX+rY+rZ,rX+1:rX+rY);
        Ezz = E2(rX+rY+1:rX+rY+rZ,rX+rY+1:rX+rY+rZ);

        %correct error
        cExx = Exx;
        cEyy = Eyy-Eyx*inv(Exx)*Exy;
        cEzz = Ezz-Ezx*inv(Exx)*Exz-(Ezy-Ezx*inv(Exx)*Exy)*inv(Eyy-Eyx*inv(Exx)*Exy)*(Eyz-Eyx*inv(Exx)*Exz);

        freqCau(m) = log(1+(trace(abs(Qxy*cEyy*conj(Qxy')))+trace(abs(Qxz*cEzz*conj(Qxz'))))/trace(abs(Qxx*cExx*conj(Qxx'))));
        %freqCau(m) = log(det(E1(1:rX,1:rX))/abs(Qxx*det(E2(1:rX,1:rX))*conj(Qxx')));

    end;

    freq = sr.*[0:cX/2-1]./cX;
end;
