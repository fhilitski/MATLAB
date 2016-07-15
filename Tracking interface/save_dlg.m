function varargout = save_dlg(varargin)
% SAVE_DLG MATLAB code for save_dlg.fig
%      SAVE_DLG, by itself, creates a new SAVE_DLG or raises the existing
%      singleton*.
%
%      H = SAVE_DLG returns the handle to a new SAVE_DLG or the handle to
%      the existing singleton*.
%
%      SAVE_DLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAVE_DLG.M with the given input arguments.
%
%      SAVE_DLG('Property','Value',...) creates a new SAVE_DLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before save_dlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to save_dlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help save_dlg

% Last Modified by GUIDE v2.5 17-Jun-2013 19:59:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @save_dlg_OpeningFcn, ...
    'gui_OutputFcn',  @save_dlg_OutputFcn, ...
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


% --- Executes just before save_dlg is made visible.
function save_dlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to save_dlg (see VARARGIN)

% Choose default command line output for save_dlg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%get input arguments
handles.Properties.traj = varargin{1};
handles.Properties.directory = varargin{2};
handles.Properties.fname = varargin{3};

[hObject,handles] = set_initial_state(hObject,handles);

guidata(hObject,handles);

% UIWAIT makes save_dlg wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function [retObject,rethandles] = set_initial_state(hObject,handles)

handles.Properties.DataFileExists = 0;
handles.Properties.RunVariableExists = 0;
handles.Properties.CalVariableExists = 0;

set(handles.edit1,'String','');
set(handles.edit2,'String','');
set(handles.edit3,'String','');

guidata(hObject,handles);
retObject = hObject;
rethandles = handles;




% --- Outputs from this function are returned to the command line.
function varargout = save_dlg_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in pushbutton1.
% --- Attemtpt to save separation
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%determine what kind of run we are saving
exp_run = get(handles.radiobutton1,'Value');
cal_run = get(handles.radiobutton2,'Value');

%if new file is necessary (no file selected)...
if ~handles.Properties.DataFileExists
    
    %get ID;
    ID = get(handles.edit1,'String');
    Day = get(handles.edit2,'String');
    
    
    FileName = [ID '.mat'];
    PathName = handles.Properties.directory;
    
    
    
else
    %if file already exists
    ID = handles.Properties.ID;
    Day = handles.Properties.Day;
    
        
    %get Filename
    PathName = handles.Properties.PathName;
    FileName = handles.Properties.FileName;
    
    
end;

%overlap lengh will always be saved as the value entered, even if its new;
Overlap_length = get(handles.edit3,'String');

%make variable name for measured separation
if exp_run
    var_name = ['R',ID,'_',Day];
end;

if cal_run
    var_name = ['R',ID,'_',Day,'cal'];
end;

traj = handles.Properties.traj;
eval([var_name ' = traj;']);

%by default, we are saving variable with name var_name into file FileName
decision = 'Yes';

%check if variable already exists in the file
if handles.Properties.DataFileExists
    try
        existing_var = load([PathName FileName], var_name);
    catch exception;
    end;
    if ~isempty(struct2cell(existing_var))
        msg_string = ['Variable ' var_name ' already exists in ' FileName '!'];
        decision = questdlg([msg_string ' Overwrite?'],'Warning!','Yes','No','Yes');
    end;
end;


if (strcmp(decision, 'Yes'))
    if handles.Properties.DataFileExists
        save([PathName FileName], var_name,'ID','Day','Overlap_length','-append','-mat');
    else
        save([PathName FileName], var_name,'ID','Day','Overlap_length','-mat');
    end;
    %Separation of two beads in is now recorded into variable with
    %consistent name
    
    set(handles.text5,'String',['Variable ' var_name ' saved into ' FileName]);
    
    if exp_run
        handles.Properties.RunVariableExists = exp_run;
    end;
    if cal_run
        
        handles.Properties.CalVariableExists = cal_run;
    end;
    handles.Properties.DataFileExists = true;
    handles.Properties.ID = ID;
    handles.Properties.Day = Day;
    handles.Properties.Overlap_length = Overlap_length;
    
    handles.Properties.PathName = PathName;
    handles.Properties.FileName = FileName;
    update_status_indicator(handles);
else
    set(handles.text5,'String','Variable has not been overwritten.');
end;
guidata(hObject,handles);


function update_status_indicator(handles)
if handles.Properties.RunVariableExists
    msg_run = 'Exists';
else
    msg_run = 'Not found';
end;
if handles.Properties.CalVariableExists
    msg_cal = 'Exists';
else
    msg_cal = 'Not found';
    
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

[FileName,PathName,FilterIndex] = uigetfile('.mat','Select MAT file with data',handles.Properties.directory);

S = {};

try
    S = load([PathName FileName], 'ID*','Day');
catch exception;
end;


if ~isempty(S)
    
    Overlap_length = 0;
    try
        vars = load([PathName FileName], 'Overlap_length');
        if ~isempty(vars) 
            Overlap_length = vars.Overlap_length;
        end;
    catch exception;
    end;
            
    
    handles.Properties.DataFileExists = true;
    
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
    if (Overlap_length ~= 0)
        set(handles.edit3,'String', Overlap_length);
    end;
    
    cal_run = {};
    data_run = {};
    
    data_run = struct2cell(load([PathName FileName], ['R',ID,'_',Day]));
    cal_run = struct2cell(load([PathName FileName], ['R',ID,'_',Day,'cal']));
    
    got_data = ~isempty(data_run);
    got_cal = ~isempty(cal_run);
    handles.Properties.RunVariableExists = got_data;
    handles.Properties.CalVariableExists = got_cal;
    
    
    update_status_indicator(handles);
    set(handles.text5,'String',['Current file: ' FileName]);
    
    
    
else
    
    
    handles.Properties.DataFileExists = 0;
    update_status_indicator(handles);
    
    
end;

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
