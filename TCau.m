function [avtimeCau,avpvalue, freq, avfreqCau, timeCau, freqCau, timepvalue] =...
    TCau(timeSeriesX,timeSeriesY,Nr,Nl,order, optwindow, sr)
%TCAU Summary of this function goes here
%   Detailed explanation goes here
%timeSeriesX Y are matrix whose every row is a variable.
%the output is the causality Y to X.
%sr is the sampling rate
data = [timeSeriesX; timeSeriesY];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate average Granger causality (avGC) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear timeCau
clear df2
fstats = 0;
for wid = 1 : size(optwindow,1)
    input_data = [];
    for Nrj = 1 : Nr
        input_data = [input_data data(:, (Nrj-1)*Nl+optwindow(wid,1):(Nrj-1)*Nl+optwindow(wid,2))];
    end  
    % causality
    [timeCau(wid),R1{wid},R2{wid},ARcoeff{wid},timepvalue(wid),covxy{wid}] = ...
        Cau(input_data(1,:), input_data(2,:), Nr, optwindow(wid,2)-optwindow(wid,1)+1, order);
    % fstats
    df1 = order;
    df2(wid) = optwindow(wid,2)-optwindow(wid,1)+1 - order - order - 1;
    fstats = fstats + (R1{wid} - R2{wid}) / R2{wid} * df2(wid) / df1;
end
avpvalue = sumF2(fstats, size(optwindow,1), df1, df2);
avtimeCau = mean(timeCau);% average GC
freq = [];
avfreqCau = [];
freqCau = [];


[rX, cX] = size(timeSeriesX);
[rY, cY] = size(timeSeriesY);

if sr ~= 0   
    lamd = 2*pi.*[0:cX/2-1]./cX;
    for wid = 1 : size(optwindow,1)
        for m = 1:length(lamd)
            clear ARmodel
            clear cov_xy
            A = reshape(ARcoeff{wid}, rX+rY, (rX+rY)*order);
            if order > 1 && size(A,3) == 1
                k = rX+rY;
                temp = A;
                clear A;
                for j = 1 : order
                    A(:,:,j) = temp(:,(j-1)*k+1:j*k);
                end
            end
            for jorder = 1 : order
                ARmodel.A(:,:,1+jorder) = A(:,:,jorder);
            end
            ARmodel.na = order;
            cov_xy = reshape(covxy{wid}, rX+rY, rX+rY);
            freqCau(wid, m) = frequencyforar(ARmodel, 1, 1, lamd(m), cov_xy, sr);
        end;
    end    
    freq = sr.*[0:cX/2-1]./cX;
    avfreqCau = mean(freqCau,1);
    freqCau = reshape(freqCau, 1, length(lamd)*size(optwindow,1));
end;

% you can also see the dynamic of the causality by timeCau and freqCau

