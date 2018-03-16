function varargout = settingoptGC(varargin)
% SETTINGOPTGC M-file for settingoptGC.fig
%      SETTINGOPTGC, by itself, creates a new SETTINGOPTGC or raises the existing
%      singleton*.
%
%      H = SETTINGOPTGC returns the handle to a new SETTINGOPTGC or the handle to
%      the existing singleton*.
%
%      SETTINGOPTGC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETTINGOPTGC.M with the given input arguments.
%
%      SETTINGOPTGC('Property','Value',...) creates a new SETTINGOPTGC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before settingoptGC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to settingoptGC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help settingoptGC

% Last Modified by GUIDE v2.5 30-Jan-2013 00:06:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @settingoptGC_OpeningFcn, ...
                   'gui_OutputFcn',  @settingoptGC_OutputFcn, ...
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


% --- Executes just before settingoptGC is made visible.
function settingoptGC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to settingoptGC (see VARARGIN)

% Choose default command line output for settingoptGC
handles.output = hObject;

% default none
set(handles.radio_wholeTS,'Value',1);
set(handles.radio_optWindow,'Value',0);
set(handles.radio_fixedwindow,'Value',0);
set(handles.fixedTW_ET,'Enable','inactive');
set(handles.fixedTW_ET,'BackgroundColor',[192/255 192/255 192/255]);
set(handles.maxTW_ET,'Enable','inactive');
set(handles.maxTW_ET,'BackgroundColor',[192/255 192/255 192/255]);
set(handles.radio_avTS,'Value',1);
set(handles.radio_spatial,'Value',0);
set(handles.index_voxels,'Enable','inactive');
set(handles.index_voxels,'BackgroundColor',[192/255 192/255 192/255]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes settingoptGC wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = settingoptGC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% Now destroy yourself
delete(hObject);

% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radio_optWindow, 'Value')
    TW_OP = 1;    % optimal window dividing
    maximumTW = str2num(get(handles.maxTW_ET, 'String'));
elseif get(handles.radio_fixedwindow,'Value')
    TW_OP = 2;     % fixed number of windows
    fixednumTW = str2num(get(handles.fixedTW_ET, 'String'));
elseif get(handles.radio_wholeTS,'Value')
    TW_OP = 3;     % none
end

if get(handles.radio_spatial, 'Value')
    SP_OP = 1;     % GCA on pairs of voxels from two ROIs, and then average the GC value
    clear Ind_Voxel
    temp = get(handles.index_voxels,'string');
    for i = 1 : size(temp,1)
        Ind_Voxel{i} = str2num(temp{i});
    end
elseif get(handles.radio_avTS,'Value')
    SP_OP = 2;     % GCA on averaged time series for two ROIs
end
handles.output = struct('TW_OP',TW_OP,'SP_OP', SP_OP);


if TW_OP == 1
    handles.output.maximumTW = maximumTW;
elseif TW_OP==2
    handles.output.fixednumTW = fixednumTW;
end

if SP_OP == 1
    handles.output.Ind_Voxel = Ind_Voxel;
end
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'We have three ways of calculating optimal Granger causality: 1) temporal--estimate GC by optimally dividing time windows; 2) spatial--estimate GC for each pair of voxels and then average the restulting GC; 3) temporal and spatial -- estimate GC for each pair of voxels by optimally dividing time window.integer number of the maximum time windows allowed to be divided by the algorithm.'},'TIPS','help','modal');
uiwait(h);


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



function fixedTW_ET_Callback(hObject, eventdata, handles)
% hObject    handle to fixedTW_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixedTW_ET as text
%        str2double(get(hObject,'String')) returns contents of fixedTW_ET as a double


% --- Executes during object creation, after setting all properties.
function fixedTW_ET_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixedTW_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'A integer number specifying the maximum number of time windows the time series is going to be optimally divided into. Recommend to be 3 or 5 for time series of moderate length. '},'TIPS','help','modal');
uiwait(h);

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'A integer number specifing the number of time windows divided by the algorithm. Recommend to be 3 or 5 for time series of moderate length. '},'TIPS','help','modal');
uiwait(h);


% --- Executes on button press in radio_wholeTS.
function radio_wholeTS_Callback(hObject, eventdata, handles)
% hObject    handle to radio_wholeTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_wholeTS
set(handles.radio_optWindow,'Value',0);
set(handles.radio_fixedwindow,'Value',0);
set(handles.fixedTW_ET,'Enable','inactive');
set(handles.fixedTW_ET,'BackgroundColor',[192/255 192/255 192/255]);
set(handles.maxTW_ET,'Enable','inactive');
set(handles.maxTW_ET,'BackgroundColor',[192/255 192/255 192/255]);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'No time window is going to be divided for GCA.'},'TIPS','help','modal');
uniwait(h);

% --- Executes on button press in radio_spatial.
function radio_spatial_Callback(hObject, eventdata, handles)
% hObject    handle to radio_spatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_spatial
set(handles.radio_avTS,'Value',0);
set(handles.index_voxels,'Enable','on');
set(handles.index_voxels,'BackgroundColor',[255/255 255/255 255/255]);

% --- Executes on button press in radio_avTS.
function radio_avTS_Callback(hObject, eventdata, handles)
% hObject    handle to radio_avTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_avTS
set(handles.radio_spatial,'Value',0);
set(handles.index_voxels,'Enable','inactive');
set(handles.index_voxels,'BackgroundColor',[192/255 192/255 192/255]);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'Granger causality on the pair of time series averaged among voxels for ROIs.'},'TIPS','help','modal');
uniwait(h);

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to index_voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of index_voxels as text
%        str2double(get(hObject,'String')) returns contents of index_voxels as a double


% --- Executes during object creation, after setting all properties.
function index_voxels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to index_voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=msgbox({'Vector specifying indeces in the data for each ROI. GCA on each pair of voxels between two ROIs'},'TIPS','help','modal');
uniwait(h);

% --- Executes during object creation, after setting all properties.
function pushbutton10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in radio_optWindow.
function radio_optWindow_Callback(hObject, eventdata, handles)
% hObject    handle to radio_optWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_optWindow
 set(handles.radio_fixedwindow,'Value',0);
 set(handles.radio_wholeTS,'Value',0);
set(handles.maxTW_ET,'Enable','on');
set(handles.maxTW_ET,'BackgroundColor',[255/255 255/255 255/255]);
set(handles.fixedTW_ET,'Enable','inactive');
set(handles.fixedTW_ET,'BackgroundColor',[192/255 192/255 192/255]);

% --- Executes on button press in radio_fixedwindow.
function radio_fixedwindow_Callback(hObject, eventdata, handles)
% hObject    handle to radio_fixedwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_fixedwindow
 set(handles.radio_optWindow,'Value',0);
 set(handles.radio_wholeTS,'Value',0);
 set(handles.fixedTW_ET,'Enable','on');
set(handles.fixedTW_ET,'BackgroundColor',[255/255 255/255 255/255]);
set(handles.maxTW_ET,'Enable','inactive');
set(handles.maxTW_ET,'BackgroundColor',[192/255 192/255 192/255]);



function index_voxels_Callback(hObject, eventdata, handles)
% hObject    handle to index_voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of index_voxels as text
%        str2double(get(hObject,'String')) returns contents of index_voxels as a double
