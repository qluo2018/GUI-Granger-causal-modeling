function varargout = GrangerCausalityGUI(varargin)
% GRANGERCAUSALITYGUI M-file for GrangerCausalityGUI.fig
%      GRANGERCAUSALITYGUI, by itself, creates a new GRANGERCAUSALITYGUI or raises the existing
%      singleton*.
%
%      H = GRANGERCAUSALITYGUI returns the handle to a new GRANGERCAUSALITYGUI or the handle to
%      the existing singleton*.
%
%      GRANGERCAUSALITYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRANGERCAUSALITYGUI.M with the given input arguments.
%
%      GRANGERCAUSALITYGUI('Property','Value',...) creates a new GRANGERCAUSALITYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GrangerCausalityGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GrangerCausalityGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GrangerCausalityGUI

% Last Modified by GUIDE v2.5 22-Aug-2013 22:31:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GrangerCausalityGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GrangerCausalityGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GrangerCausalityGUI is made visible.
function GrangerCausalityGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GrangerCausalityGUI (see VARARGIN)

% Choose default command line output for GrangerCausalityGUI
global TimeS;
global totalLength;

TimeS = [];
totalLength = 0;
handles.output = hObject;
handles.GC_OP = 1;
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes GrangerCausalityGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GrangerCausalityGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
pathDisplay = pwd;
% n = length(pathDisplay);
% pathDisplay = pathDisplay(n-26:n);
set(handles.path_T,'String',pathDisplay);



%% TEXT callback and create function
function path_T_Callback(hObject, eventdata, handles)
% hObject    handle to path_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path_T as text
%        str2double(get(hObject,'String')) returns contents of path_T as a double



% --- Executes during object creation, after setting all properties.
function path_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

function Nv_ET_Callback(hObject, eventdata, handles)
% hObject    handle to Nv_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nv_ET as text
%        str2double(get(hObject,'String')) returns contents of Nv_ET as a double


% --- Executes during object creation, after setting all properties.
function Nv_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nv_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ns_ET_Callback(hObject, eventdata, handles)
% hObject    handle to Ns_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ns_ET as text
%        str2double(get(hObject,'String')) returns contents of Ns_ET as a double


% --- Executes during object creation, after setting all properties.
function Ns_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ns_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SR_ET_Callback(hObject, eventdata, handles)
% hObject    handle to SR_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SR_ET as text
%        str2double(get(hObject,'String')) returns contents of SR_ET as a double


% --- Executes during object creation, after setting all properties.
function SR_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SR_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Nr_ET_Callback(hObject, eventdata, handles)
% hObject    handle to Nr_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nr_ET as text
%        str2double(get(hObject,'String')) returns contents of Nr_ET as a double
global totalLength;

Nr_str = get(hObject,'String');
Ns = totalLength/str2double(Nr_str);
set(handles.Ns_ET,'String',Ns);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Nr_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nr_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bandEdge_ET_Callback(hObject, eventdata, handles)
% hObject    handle to bandEdge_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bandEdge_ET as text
%        str2double(get(hObject,'String')) returns contents of bandEdge_ET as a double


% --- Executes during object creation, after setting all properties.
function bandEdge_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bandEdge_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function downSampling_ET_Callback(hObject, eventdata, handles)
% hObject    handle to downSampling_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of downSampling_ET as text
%        str2double(get(hObject,'String')) returns contents of downSampling_ET as a double


% --- Executes during object creation, after setting all properties.
function downSampling_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to downSampling_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upSampling_ET_Callback(hObject, eventdata, handles)
% hObject    handle to upSampling_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upSampling_ET as text
%        str2double(get(hObject,'String')) returns contents of upSampling_ET as a double


% --- Executes during object creation, after setting all properties.
function upSampling_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upSampling_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numDif_ET_Callback(hObject, eventdata, handles)
% hObject    handle to numDif_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numDif_ET as text
%        str2double(get(hObject,'String')) returns contents of numDif_ET as a double


% --- Executes during object creation, after setting all properties.
function numDif_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numDif_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function modelOrder_ET_Callback(hObject, eventdata, handles)
% hObject    handle to modelOrder_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modelOrder_ET as text
%        str2double(get(hObject,'String')) returns contents of modelOrder_ET as a double


% --- Executes during object creation, after setting all properties.
function modelOrder_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelOrder_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numS_ET_Callback(hObject, eventdata, handles)
% hObject    handle to numS_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numS_ET as text
%        str2double(get(hObject,'String')) returns contents of numS_ET as a double


% --- Executes during object creation, after setting all properties.
function numS_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numS_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Buttons CallBack functions

% --- Executes on button press in open_B.
function open_B_Callback(hObject, eventdata, handles)
% hObject    handle to open_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global TimeS;
global totalLength;
clc;

[filename,pathname] = uigetfile;

if (filename == 0 & pathname == 0)
    msgbox('Please Load Data!','Data Error');
    return;
else
    TimeSS = load([pathname, filename]);
    FNs = fieldnames(TimeSS);
    TimeS = getfield(TimeSS,FNs{1});
    [Nv,Ns] = size(TimeS);
    totalLength = Ns;
    for i=1:min(Nv,5)
        subplot(min(Nv,5),1,i); plot(TimeS(i,:));hold on;       
    end
    wholePath = [pathname,filename];
    n=length(wholePath);
    pathDisplay = wholePath(n-26:n);
    set(handles.path_T,'String',pathDisplay);
    set(handles.Nv_ET,'String',num2str(Nv));
    set(handles.Ns_ET,'String',num2str(Ns));
    set(handles.Nr_ET,'String',num2str(1));
    set(handles.SR_ET,'String',num2str(1000));
    set(handles.downSampling_ET,'String',num2str(1000));
    set(handles.upSampling_ET,'String',num2str(1000));
    %set(handles.dataLength,'String',num2str(Ns));
end
guidata(hObject, handles);

% --- Executes on button press in DItip_B.
function DItip_B_Callback(hObject, eventdata, handles)
% hObject    handle to DItip_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'Nv = Number of Variables'; '';'Ns = Number of time points'; '';'Nr = Number of realizations';'';'SR = sampling rate (Hz, =1/TR)'},'TIPS','help','modal');
uiwait(h);

% --- Executes on button press in Ftips_B.
function Ftips_B_Callback(hObject, eventdata, handles)
% hObject    handle to Ftips_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=msgbox({'A vector of frequency points, in ascending order, in the range 0 to 1. The value 1 corresponds to half the sample frequency. This vector must have even length. Tunable.'},'TIPS','help','modal');
uiwait(h);


% --- Executes on button press in SRMTips_B.
function SRMTips_B_Callback(hObject, eventdata, handles)
% hObject    handle to SRMTips_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'Enter the sampling rate after modification in Hz. If this number is equal to the sampling rate of the data, nothing happens'},'TIPS','help','modal');
uiwait(h);


% --- Executes on button press in ARIMATips_B.
function ARIMATips_B_Callback(hObject, eventdata, handles)
% hObject    handle to ARIMATips_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'Set the order of ARIMA model. The ARIMA model will be same as AR model by setting the number of differences to 0. The BEKK model order specifies the model order in the variance in GCSDN'},'TIPS','help','modal');
uiwait(h);


% --- Executes on button press in start_B.
function start_B_Callback(hObject, eventdata, handles)
% hObject    handle to start_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global TimeS;
global totalLength;

Nv = str2num(get(handles.Nv_ET,'String'));
Ns = str2num(get(handles.Ns_ET,'String'));
Nr = str2num(get(handles.Nr_ET,'String'));
SR = str2num(get(handles.SR_ET,'String'));
BE = str2num(get(handles.bandEdge_ET,'string'));
sampling_OP = get(handles.downSampling_RB,'Value');
GC_OP = handles.GC_OP;

dsF = str2num(get(handles.downSampling_ET,'String'));
usF = str2num(get(handles.upSampling_ET,'String'));
NoD = str2num(get(handles.numDif_ET,'String'));
order = str2num(get(handles.modelOrder_ET,'String'));
bekk_order = str2num(get(handles.bekkmodel_order_ET,'String'));
N_sim = str2num(get(handles.numS_ET,'String'));
totalLoops = N_sim;
progress = 0;
if ~isempty(BE)
    totalLoops=totalLoops+length(BE)/2*Nv;
end;
if sampling_OP == 1
    if dsF ~= SR
        totalLoops=totalLoops+Nr;
    end;
    fre = dsF;
else
    if usF ~= SR
        totalLoops=totalLoops+Nr;
    end;
    fre = upF;
end;
if NoD > 0
    totalLoops = totalLoops+Nr;
end;
totalLoops = totalLoops+Nv*Nv;

progressbar;
%% filter
if ~isempty(BE)
    for i=1:length(BE)/2
        [w1,w2] = cheby1(3,0.5,BE(2*(i-1)+1:2*i),'stop');
        for j=1:Nv
            timeSf(j,:) = filter(w1,w2,TimeS(j,:));
            progress = progress+1;
            progressbar(progress/totalLoops);
        end;
    end;
else
    timeSf = TimeS;
end;

%% sampling
dtimeS = [];
if sampling_OP == 1
    dsRate = floor(SR/dsF);
    if dsRate ~= 1
        for i=1:Nr
            for j=1:floor(Ns/dsRate)
                dtimeS(:,j+(i-1)*floor(Ns/dsRate)) = timeSf(:,j*dsRate+(i-1)*Ns);
            end;
            progress = progress+1;
            progressbar(progress/totalLoops);
        end;
    else
        dtimeS = timeSf;
    end;
else
    dsRate = SR/usF;
    if dsRate ~= 1;
        xx = [1:dsRate:Ns];
        for i=1:Nr
            dtimeS = [dtimeS,spline(1:Ns,timeSf(:,1+Ns*(i-1):Ns*i),xx)];
            progress = progress+1;
            progressbar(progress/totalLoops);
        end;
    else
        dtimeS = timeSf;
    end;
end;
Ns = Ns/dsRate;

%% ARIMA
dtimeSA = [];
if (NoD > 0)
    for i=1:Nr
        dtimeSA = [dtimeSA,integrateDifference(dtimeS(:,1+(i-1)*Ns:i*Ns),NoD)];
        progress = progress+1;
        progressbar(progress/totalLoops);
    end;      
else
    dtimeSA = dtimeS;
end;

Ns = Ns-NoD*Nr;

%% boostrapping
if GC_OP == 3
    totalLoops = Nv*Nv * Nv;
    TGC_setting = handles.TGC_setting;
    if TGC_setting.TW_OP == 1
        maximumTW = TGC_setting.maximumTW;
    elseif TGC_setting.TW_OP == 2
        fixednumTW = TGC_setting.fixednumTW;
    end
    if TGC_setting.SP_OP == 1
        Ind_Voxel = TGC_setting.Ind_Voxel;
    end
    causMatrix = zeros(Nv,Nv);
    pvalMatrix = eye(Nv);
    clear timeCau;
    for j=1:Nv
        for k=1:Nv
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % estimating causality for each direction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if j ~= k               
                withinROI = 0;
                if TGC_setting.SP_OP == 1                    
                    % if SP_OP option is on, then we can ignore some
                    % pairs of voxels within the same ROI
                    for roiind = 1 : size(Ind_Voxel,2)
                        if ismember(j, Ind_Voxel{roiind}) && ismember(k,  Ind_Voxel{roiind}) 
                            withinROI = 1;
                            break
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % optimally dividing time windows
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if withinROI == 0
                    if k > j
                        % temporal
                        if TGC_setting.TW_OP == 1
                            % optimal time window for each pair of TS in the
                            % data set
                            clear optchangepoint
                            [optcoeff, opterror, optchangepoint, BIC, optexit] = ...
                                opttvCau(dtimeSA([j,k],:),Nr,Ns,order, 10, maximumTW);
                            %optchangepoint
                            for i = 1 : size(optchangepoint, 1)
                                optwindow{j,k}(i,:) = [optchangepoint(i,1),optchangepoint(i,2)];
                            end
                            if size(optchangepoint, 1) == 0
                                optwindow{j,k}(1,:) = [1,Ns];
                            end
                        elseif TGC_setting.TW_OP == 2
                            % fixed time window for each pair of TS in the data set
                            clear optchangepoint
                            windowlength = floor( (Ns-1) / fixednumTW);
                            if windowlength > 4
                                optwindow{j,k}(1,:) = [1, windowlength];
                                for i = 2 : fixednumTW-1
                                    optwindow{j,k}(i,:) =  [windowlength*(i-1)+1, windowlength*i];
                                end
                                optwindow{j,k}(fixednumTW,:) =  [windowlength*(fixednumTW-1)+1,Ns];
                            else
                                disp('window is too short to be defined. It has been ignored.')
                                optwindow{j,k}(1,:) = [1,Ns];
                            end
                        else
                            optwindow{j,k}(1,:) = [1,Ns];
                        end
                    else
                        optwindow{j,k} = optwindow{k,j};
                    end                
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % causality estimation from kth region to jth region
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [TcausMatrix(j,k), TpvalMatrix(j,k), f, fc, TtimeCau{j,k}, timefc, Ttimepvalue{j,k}] = ...
                        TCau(dtimeSA(j,:),dtimeSA(k,:),Nr,Ns,order,optwindow{j,k}, 0);
                end
            end
            progress = progress+1;
            progressbar(progress/totalLoops);
        end
    end      
    clear causMatrix
    clear stdMatrix
    clear pvalMatrix
    clear xLables
    % spatial
    if TGC_setting.SP_OP == 1
        % voxel
        figure;
        count = 1;
        numROI = size(Ind_Voxel,2);
        for j = 1 : numROI
            for k = 1 : numROI
                if j ~= k
                    clear indexX
                    clear indexY
                    indexJ = Ind_Voxel{j}; 
                    indexK= Ind_Voxel{k};
                    % from k to j
                    causMatrix(count) = mean(mean(TcausMatrix(indexJ,indexK)));
                    clear temp
                    temp = TcausMatrix(indexJ,indexK);
                    stdMatrix(count) = std(temp(:));
                    % display                    
                    xLables{count} = [num2str(j),'-->',num2str(k)];            
                    count = count + 1;
                end
            end
        end 
        errorbar(1:count-1, causMatrix, stdMatrix, 'r-x');
        xticklabel_rotate(1:count-1,70,xLables);
        ylabel('causality')
        xlabel('ROI A --> ROI B')   
        colheads = ['Direction           ';'   mean causality   ';'      std causality '];
        disCauTable([causMatrix; stdMatrix]',colheads, xLables, TGC_setting, GC_OP)
    else
        figure;        
        % none 
        count1 = 0;
        count = 0;
        for j=1:Nv
            for k=1:Nv
                count1 = count1 + 1;
                if j ~= k
                    count = count + 1;
                    causMatrix(count) = TcausMatrix(j,k);
                    pvalMatrix(count) = TpvalMatrix(j,k);
                    xLables{count} = [num2str(j),'-->',num2str(k)];      
                    subplot(Nv, Nv, count1);
                    plot(TtimeCau{j,k},'b-+');
                    hold on
                    plot(Ttimepvalue{j,k}, 'r-x')
                    hold off
                    if count == 1
                        ylabel('Granger Causality')
                        legend('Granger causality', 'p-value', 'Location','best')                        
                    end
                    title([num2str(k), '-->', num2str(j)])
%                     legend('casality','pvalue')
                    xlim([0.5, size(TtimeCau{j,k},2)+0.5]);
                end
            end
        end
        figure
        subplot(1,2,1)
        bar(causMatrix)
        ylabel('causality')
        xticklabel_rotate(1:count,70, xLables);
        subplot(1,2,2)
        bar(pvalMatrix)
        ylabel('p-value')
        xticklabel_rotate(1:count,70, xLables);
        xlabel('ROI A --> ROI B') 
        colheads = ['Direction           ';'      causality     ';'             p-value'];
        disCauTable([causMatrix; pvalMatrix]', colheads, xLables, TGC_setting, GC_OP)
    end
elseif GC_OP == 4
    if Nr > 1
        % SDN identification
        figure('Name','SDN identification')
        [cc0,pp0] = BOLD_SDN_Identify(dtimeSA, Nr, Ns, order);
    end    
    
    % GCSDN
    progressbar
    progress = 0;
    totalLoops = Nv * (Nv-1) / 2;
    clear timefreqCau
    fd.EDFreq = fre/2; fd.STFreq = 0;
    fd.NFFT = 256;  fd.fs = fre;
    for j = 1 : Nv-1
        for k = j+1 : Nv
            clear input_data
            input_data.timeseriesdata = dtimeSA([j,k],:)';
            input_data.Nr = Nr;
            input_data.Nl = Ns;
            index_X = 1; % it can be multivariate
            index_Y = 2; % it can be multivariate
            outputarmabekk = mv_grangerarmabekk4Repeat(input_data, ...
                order, bekk_order, index_X, index_Y, fd);
            TcausMatrix(k,j) = outputarmabekk.granger(1);  % j to k % Log-likelihood ratio version
            TpvalMatrix(k,j) = outputarmabekk.granger(3);  % j to k
            VcausMatrix(k,j) = outputarmabekk.granger(5);  % j to k  % variance ratio version
            
            TcausMatrix(j,k) = outputarmabekk.granger(2);  % k to j
            TpvalMatrix(j,k) = outputarmabekk.granger(4);  % k to j
            VcausMatrix(j,k) = outputarmabekk.granger(6); % k to j
            
            timefreqCau{k,j} = outputarmabekk.fgranger;   % frequency domain result
            
            progress = progress+1;
            progressbar(progress/totalLoops);            
        end        
    end
    
    
    count = 1;
    for j=1:Nv
        for k=1:Nv
            if j ~= k
                causMatrix(count) = TcausMatrix(j,k);
                pvalMatrix(count) = TpvalMatrix(j,k);                                
                xLables{count} = [num2str(j),'-->',num2str(k)];
                count = count + 1;
            end;
        end;
    end;
    
    figure
    subplot(1,2,1)
    bar(causMatrix)
    ylabel('causality')
    xticklabel_rotate(1:count-1,70, xLables);
    subplot(1,2,2)
    bar(pvalMatrix)
    ylabel('p-value')
    xticklabel_rotate(1:count-1,70, xLables);
    xlabel('ROI A --> ROI B')
    colheads = ['Direction           ';'      causality     ';'             p-value'];
    TGC_setting = 0;    
    disCauTable([causMatrix; pvalMatrix]', colheads, xLables, TGC_setting, GC_OP)

else
    [A,E] = armorf(dtimeSA,Nr,Ns,order);
    for simC = 1:N_sim        
        gData = samplingAR(A,E,dtimeSA,Nv,order,Nr,Ns);        
        for j=1:Nv
            for k=1:Nv
                if j ~=k
                    cData = [1:Nv];
                    cData([j,k]) = [];
                    if GC_OP == 1;
                        causMatrix{j,k,simC}= CCau(gData(j,:),gData(k,:),gData(cData,:),Nr,Ns,order,0);
                    elseif GC_OP == 2
                        causMatrix{j,k,simC}= PCau(gData(j,:),gData(k,:),gData(cData,:),Nr,Ns,order,0);                   
                    end;
                end;
            end;
        end;        
        progress = progress+1;
        progressbar(progress/totalLoops);        
    end;  
    
    xLables{1} = '';
    count = 1;
    for j=1:Nv
        for k=1:Nv
            if j ~=k
                meanArray(count) = mean(cell2mat(causMatrix(j,k,:)));
                varArray(count) = sqrt(var(cell2mat(causMatrix(j,k,:))));
                count = count + 1;
                xLables{count} = [num2str(k),'-->',num2str(j)];
            end;
        end;
    end;
    
    figure();
    errorbar([1:nchoosek(Nv,2)*2],meanArray,1.96*varArray);
    set(gca,'xtick',[0:count-1]);
    set(gca,'xticklabel',xLables);
    title('Granger Causality in time domain');
    ylabel('Granger Causality');
    %     xticklabel_rotate([],70);
    colheads = ['Direction           ';'   mean causality   ';'      std causality '];
    TGC_setting = 0;    
    for xi = 2 : length(xLables)        
        xLables_temp{xi-1} = xLables{xi};
    end  
    disCauTable([meanArray; varArray]', colheads, xLables_temp, TGC_setting, GC_OP)
end


%% Frequency
deData = dtimeSA;
  
if GC_OP ~= 3
    figure;
end

% ylimt = 0;
% for j=1:Nv
%     if GC_OP ~= 3       
%         for i=1:Nr
%             pwelch(deData(j,1+Ns*(i-1):Ns*i),[],[],[],fre);
%             [pxx(:,i),f] = pyulear(deData(j,1+Ns*(i-1):Ns*i),order,[],fre);
%         end;
%         if max(max(pxx)) > ylimt
%             ylimt = max(max(pxx));
%         end
%     end
% end
for j=1:Nv
    if GC_OP ~= 3       
        subplot(Nv,Nv,j+Nv*(j-1));
        for i=1:Nr
            pwelch(deData(j,1+Ns*(i-1):Ns*i),[],[],[],fre);
            %[pxx(:,i),f] = pyulear(deData(j,1+Ns*(i-1):Ns*i),order,[],fre);
        end;
        title([num2str(j)]);
        %plot(f, mean(pxx,2));
        %ylim([0, ylimt])
    end    
    %plot(f*60,fc(1:end/2));
    %axis([0 0.026 0 0.4]);    
    for k=1:Nv
        if j ~=k
            cData = [1:Nv];
            cData([j,k]) = [];
            
            if GC_OP == 1
                [c,f,fc] = CCau(deData(j,:),deData(k,:),deData(cData,:),Nr,Ns,order,fre);
            elseif GC_OP == 2
                [c,f,fc] = PCau(deData(j,:),deData(k,:),deData(cData,:),Nr,Ns,order,fre);
            elseif GC_OP == 3 && TGC_setting.SP_OP == 0
                [c,p, f,fc, tc, timefreqCau{j,k}] = TCau(deData(j,:),deData(k,:),Nr,Ns,order,optwindow{j,k},fre);              
            end 
            if GC_OP ==1 || GC_OP == 2
                subplot(Nv,Nv,(j-1)*Nv+k);
                plot(f,fc);
                %ylim([0 ylimt])
                
                %axis([0 fre/2 0 0.21]);
                title([num2str(k),'-->',num2str(j)]);
            elseif GC_OP == 4
                subplot(Nv,Nv,(j-1)*Nv+k);    
                if j < k          
                    f = timefreqCau{k,j}(1,:);
                    fc = timefreqCau{k,j}(2,:);                                      
                else
                    f = timefreqCau{j,k}(1,:);
                    fc = timefreqCau{j,k}(3,:); 
                end
                plot(f,fc); 
                %ylim([0 ylimt])
                %axis([0 fre/2 0 0.21]);
                title([num2str(k),'-->',num2str(j)]);
                count =  count + 1;
            end
        end;

    end;
end;

if GC_OP == 3 && TGC_setting.SP_OP == 0
    figure();
    count = 0;
    for j=1:Nv
        for k = 1 : Nv
            if k ~= j
                count = count + 1;
                subplot(Nv, Nv, count);
                x_time = [];
                y_freq = [];
                z_freqCau = [];
                z_freqCau = timefreqCau{j,k};
                numWindow = size(optwindow{j,k},1);
                numFreq = size(f,2);
                tt = 1;
                for fi = 1 : numFreq
                    for ti = 1 : numWindow
                        x_time(tt) = ti;
                        y_freq(tt) = f(fi);
                        tt = tt + 1;
                    end
                end
                plot3(x_time, y_freq, z_freqCau, 'b-')
                title([num2str(j),'-->',num2str(k)]);
            end
        end
    end
    xlabel('Time')
    ylabel('Frequency')
    zlabel('Causality')
end 
progressbar(1);

save('myresult.mat')

%% Radio Buttons

% --- Executes on button press in downSampling_RB.
function downSampling_RB_Callback(hObject, eventdata, handles)
% hObject    handle to downSampling_RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of downSampling_RB
set(hObject,'Value',1);
set(handles.upSampling_RB,'Value',0);
set(handles.upSampling_ET,'Enable','off');
set(handles.upSampling_ET,'BackgroundColor',[192/255 192/255 192/255]);
set(handles.downSampling_ET,'Enable','on');
set(handles.downSampling_ET,'BackgroundColor',[255/255 255/255 255/255]);

% --- Executes on button press in upSampling_RB.
function upSampling_RB_Callback(hObject, eventdata, handles)
% hObject    handle to upSampling_RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of upSampling_RB
set(hObject,'Value',1);
set(handles.downSampling_RB,'Value',0);
set(handles.downSampling_ET,'Enable','off');
set(handles.downSampling_ET,'BackgroundColor',[192/255 192/255 192/255]);
set(handles.upSampling_ET,'Enable','on');
set(handles.upSampling_ET,'BackgroundColor',[255/255 255/255 255/255]);

function maxTW_ET_Callback(hObject, eventdata, handles)
% hObject    handle to maxTW_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxTW_ET as text
%        str2double(get(hObject,'String')) returns contents of maxTW_ET as a double


% --- Executes during object creation, after setting all properties.
function maxTW_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxTW_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'A integer number of the maximum time windows allowed to be divided by the algorithm. Recommend to be 3 or 5 for time series of moderate length. '},'TIPS','help','modal');
uiwait(h);


% --- Executes when selected object is changed in uipanel12.
function uipanel12_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel12 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.CGC_RB,'Value')
    handles.GC_OP = 1;
    set(handles.bekkmodel_order_ET,'Enable','off');
    set(handles.bekkmodel_order_ET,'BackgroundColor',[192/255 192/255 192/255]);
    set(handles.numS_ET,'Enable','on');
    set(handles.numS_ET,'BackgroundColor','w');
elseif get(handles.PGC_RB,'Value')
    handles.GC_OP = 2;
    set(handles.bekkmodel_order_ET,'Enable','off');
    set(handles.bekkmodel_order_ET,'BackgroundColor',[192/255 192/255 192/255]);
    set(handles.numS_ET,'Enable','on');
    set(handles.numS_ET,'BackgroundColor','w');
elseif get(handles.TGC_RB,'Value')
    handles.GC_OP = 3;
    set(handles.numS_ET,'Enable','off');
    set(handles.numS_ET,'BackgroundColor',[192/255 192/255 192/255]);
    tempoutput = settingoptGC; 
    handles.TGC_setting = tempoutput;
    set(handles.bekkmodel_order_ET,'Enable','off');
    set(handles.bekkmodel_order_ET,'BackgroundColor',[192/255 192/255 192/255]);
elseif get(handles.GCSDN_RB,'Value')
    handles.GC_OP = 4;
    set(handles.numS_ET,'Enable','off');
    set(handles.numS_ET,'BackgroundColor',[192/255 192/255 192/255]);
    set(handles.bekkmodel_order_ET,'Enable', 'on');
    set(handles.bekkmodel_order_ET,'BackgroundColor','w');    
    
end
guidata(hObject, handles);




% --- Executes on key press with focus on DItip_B and none of its controls.
function DItip_B_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to DItip_B (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'Data matrix with each row for one variable each column for one time point'},'TIPS','help','modal');
uiwait(h);

% --- Executes on key press with focus on pushbutton10 and none of its controls.
function pushbutton10_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function bekkmodel_order_ET_Callback(hObject, eventdata, handles)
% hObject    handle to bekkmodel_order_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bekkmodel_order_ET as text
%        str2double(get(hObject,'String')) returns contents of bekkmodel_order_ET as a double


% --- Executes during object creation, after setting all properties.
function bekkmodel_order_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bekkmodel_order_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
