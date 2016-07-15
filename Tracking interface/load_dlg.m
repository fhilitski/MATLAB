function varargout = load_dlg(varargin)
% LOAD_DLG MATLAB code for load_dlg.fig
%      LOAD_DLG, by itself, creates a new LOAD_DLG or raises the existing
%      singleton*.
%
%      H = LOAD_DLG returns the handle to a new LOAD_DLG or the handle to
%      the existing singleton*.
%
%      LOAD_DLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_DLG.M with the given input arguments.
%
%      LOAD_DLG('Property','Value',...) creates a new LOAD_DLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before load_dlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to load_dlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help load_dlg

% Last Modified by GUIDE v2.5 19-Mar-2014 16:57:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @load_dlg_OpeningFcn, ...
    'gui_OutputFcn',  @load_dlg_OutputFcn, ...
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


% --- Executes just before load_dlg is made visible.
function load_dlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to load_dlg (see VARARGIN)

% Choose default command line output for load_dlg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.Properties.directory = varargin{2};

[hObject,handles] = set_initial_state(hObject,handles);
guidata(hObject,handles);

% UIWAIT makes load_dlg wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function [retObject,rethandles] = set_initial_state(hObject,handles)

handles.Properties.DataFileExists = 0;
handles.Properties.RunVariableExists = 0;
handles.Properties.CalVariableExists = 0;
handles.Properties.Plots = 0;

set(handles.edit1,'String','');
set(handles.edit2,'String','');
set(handles.edit3,'String','');
set(handles.pushbutton4,'Enable','off');
set(handles.pushbutton5,'Enable','off');
set(handles.pushbutton6,'Enable','off');
set(handles.pushbutton7,'Enable','off');
set(handles.pushbutton8,'Enable','off');

guidata(hObject,handles);
retObject = hObject;
rethandles = handles;




% --- Outputs from this function are returned to the command line.
function varargout = load_dlg_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2



function update_status_indicator(handles)

%enable or disble plot buttons based on separation availability
if handles.Properties.RunVariableExists
    msg_run = 'Exists';
    set(handles.pushbutton5,'Enable','on');
else
    msg_run = 'Not found';
    set(handles.pushbutton5,'Enable','off');
end;

if handles.Properties.CalVariableExists
    msg_cal = 'Exists';
    set(handles.pushbutton6,'Enable','on');
else
    msg_cal = 'Not found';
    set(handles.pushbutton6,'Enable','off');
end;

%enable or disble FT button based on any separation availability
if handles.Properties.RunVariableExists || handles.Properties.CalVariableExists
    set(handles.pushbutton7,'Enable','on');
    set(handles.pushbutton8,'Enable','on');
else
    set(handles.pushbutton7,'Enable','off');
    set(handles.pushbutton8,'Enable','off');
end;

set(handles.text8, 'String', msg_run);
set(handles.text9, 'String', msg_cal);

%check if both calibration and run are available,
%and update status of Ubmrella button
if (handles.Properties.RunVariableExists && handles.Properties.CalVariableExists)
    set(handles.pushbutton4,'Enable','on');
else
    set(handles.pushbutton4,'Enable','off');
end;


% --- Executes on button press in pushbutton3.
% Select .mat file with data

function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[hObject, handles] = set_initial_state(hObject,handles);

[FileName,PathName,FilterIndex] = uigetfile('.mat','Select MAT file with data',...
    handles.Properties.directory);

S = {};

try
    S = load([PathName FileName], 'ID','Day');
catch exception;
end;


if ~isempty(S)
    
    handles.Properties.DataFileExists = true;
    Overlap_length = 0;
    
    try
        vars = load([PathName FileName], 'Overlap_length');
        if ~isempty(vars)
            Overlap_length = vars.Overlap_length;
        end;
    catch exception;
    end;
    
    ID = S.ID;
    %old versions of software had ID with underscore (ex. run05_).
    %newer versions do not require underscore (ex. run05),
    %and add underscore automatically.
    %check if ID is on the old format and make it new
    if strcmp(ID(end), '_')
        ID = ID(1:end-1);
    end;
    
    Day = S.Day;
    handles.Properties.Day = Day;
    handles.Properties.ID = ID;
    handles.Properties.PathName = PathName;
    handles.Properties.FileName = FileName;
    handles.Properties.Overlap_length = Overlap_length;
    
    
    set(handles.edit1,'String', ID);
    set(handles.edit2,'String', Day);
    set(handles.edit3,'String', Overlap_length);
    
    cal_run = {};
    data_run = {};
    
    
    data_run = struct2cell(load([PathName FileName], ['R',ID,'_',Day]));
    cal_run = struct2cell(load([PathName FileName], ['R',ID,'_',Day,'cal']));
    
    got_data = ~isempty(data_run);
    got_cal = ~isempty(cal_run);
    handles.Properties.RunVariableExists = got_data;
    handles.Properties.CalVariableExists = got_cal;
    
    
    update_status_indicator(handles);
    
else
    
    handles.Properties.DataFileExists = 0;
    update_status_indicator(handles);
    
    
end;

set(handles.text5,'String',['Current file: ' FileName]);
guidata(hObject,handles);


% --- Executes on button press in pushbutton4.
% --- Perform Ubmrella sampling
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text5,'String','Performing Umberlla sampling!');

%generate filenames
ID = handles.Properties.ID;
Day = handles.Properties.Day;
PathName = handles.Properties.PathName;
FileName = handles.Properties.FileName;

exp_name = ['R',ID,'_',Day];
cal_name = ['R',ID,'_',Day,'cal'];
u_name = ['U',ID,'_',Day];
ru_name = ['Ru',ID,'_',Day];


load([PathName FileName]);

eval(['experiment = ' exp_name ';']);
eval(['calibration = ' cal_name ';']);

[n,h] = get_optimal_binning(experiment);

command = ['[' u_name ',' ru_name '] = Umbrella_new(' exp_name ',' cal_name ',' num2str(h) ');'];
%eval(['[U',ID,Day,',Ru',ID,Day,'] = Umbrella(' experiment ',' calibration ',' h ');']);
eval(command);
eval(['u = ' u_name ';']);
eval(['ru = ' ru_name ';']);

if (( min(min(u)) ~= 0 ) && ( min(min(ru)) ~= 0))
    %append results to mat file
    save([PathName FileName], u_name, ru_name ,'-append','-mat');
    set(handles.text5,'String',['Results are saved into ' PathName FileName '!']);
else
    set(handles.text5,'String','Umbrella was not performed!');
end;


% --- Executes on button press in pushbutton5.
% --- Plot Separation in separate window;
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check if there is already a plot of calibration or spearation;
if (handles.Properties.Plots == 0)
    hPlots = figure;
    handles.Properties.Plots = hPlots;
else
    hPlots = handles.Properties.Plots;
end;


%generate filenames
ID = handles.Properties.ID;
Day = handles.Properties.Day;
PathName = handles.Properties.PathName;
FileName = handles.Properties.FileName;

exp_name = ['R',ID,'_',Day];

load([PathName FileName]);

eval(['experiment = ' exp_name ';']);


figure(hPlots);
hold on;
h_plot_exp = plot(experiment);
title_msg = ['Separation between beads: experiment.'];
axis tight;
set(h_plot_exp,'LineStyle','none','Marker','.','Color','red');
title(title_msg);
xlabel('Frame #');
ylabel('Distance, \mu m');
hold off;

%update guidata
guidata(hObject, handles);



% --- Executes on button press in pushbutton6.
% --- Plot calibration in the same window as experimental plot;
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%generate filenames

if (handles.Properties.Plots == 0)
    hPlots = figure;
    handles.Properties.Plots = hPlots;
else
    hPlots = handles.Properties.Plots;
end;


ID = handles.Properties.ID;
Day = handles.Properties.Day;
PathName = handles.Properties.PathName;
FileName = handles.Properties.FileName;

cal_name = ['R',ID,'_',Day,'cal'];

load([PathName FileName]);

eval(['calibration = ' cal_name ';']);

figure(hPlots);
hold on;
h_plot_exp = plot(calibration);
title_msg = ['Separation between beads: experiment+calibration.'];
axis tight;
set(h_plot_exp,'LineStyle','none','Marker','.','Color','blue');
title(title_msg);
xlabel('Frame #');
ylabel('Distance, \mu m');
hold off;

guidata(hObject,handles);




function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when FT button is pressed;
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%generate filenames
ID = handles.Properties.ID;
Day = handles.Properties.Day;
PathName = handles.Properties.PathName;
FileName = handles.Properties.FileName;

exp_name = ['R',ID,'_',Day];
cal_name = ['R',ID,'_',Day,'cal'];

load([PathName FileName]);

%variable named experiment has experimental separation
if handles.Properties.RunVariableExists
    eval(['experiment = ' exp_name ';']);
    fc_exp = fourier_transform_routine(experiment);
    
end;

%variable named calibration has calibration separation
if handles.Properties.CalVariableExists
    eval(['calibration = ' cal_name ';']);
    fc_cal = fourier_transform_routine(calibration);
end;


% --- Executes on button press Separation MSD
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%generate filenames
ID = handles.Properties.ID;
Day = handles.Properties.Day;
PathName = handles.Properties.PathName;
FileName = handles.Properties.FileName;

exp_name = ['R',ID,'_',Day];
cal_name = ['R',ID,'_',Day,'cal'];

load([PathName FileName]);

%variable named experiment has experimental separation
if handles.Properties.RunVariableExists
    eval(['experiment = ' exp_name ';']);
    fc_exp = msd_plot_routine(experiment);
    
end;

%variable named calibration has calibration separation
if handles.Properties.CalVariableExists
    eval(['calibration = ' cal_name ';']);
    fc_cal = msd_plot_routine(calibration);
end;
