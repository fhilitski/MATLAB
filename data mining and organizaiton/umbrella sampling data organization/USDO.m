function varargout = USDO(varargin)
% USDO MATLAB code for USDO.fig
%      ************************************     
%      |                                  |                  
%      | Umbrella Sampling Data Organizer |
%      |                                  | 
%      ************************************

%      USDO, by itself, creates a new USDO or raises the existing
%      singleton*.
%
%      H = USDO returns the handle to a new USDO or the handle to
%      the existing singleton*.
%
%      USDO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in USDO.M with the given input arguments.
%
%      USDO('Property','Value',...) creates a new USDO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before USDO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to USDO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help USDO

% Last Modified by GUIDE v2.5 21-Jan-2015 11:37:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @USDO_OpeningFcn, ...
    'gui_OutputFcn',  @USDO_OutputFcn, ...
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

%Umbrella Sampling Data Organization (USDO) application


% --- Executes just before USDO is made visible.
function USDO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to USDO (see VARARGIN)

% Choose default command line output for USDO
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set_initial_state(hObject, handles);

% UIWAIT makes USDO wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function set_initial_state(hObject, handles)
%Restore initial state of the gui

% Define necessary variables
% Current working folder
handles.Current_Folder = '';
handles.Data_List = {};
handles.All_Figures = {};
set(handles.text2,'String',handles.Current_Folder);
set(handles.uipanel1,'Visible','off');
set(handles.listbox1,'String',{});
set(handles.listbox1,'Value',1);
set(handles.checkbox1,'Value',0);

drawnow;

guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = USDO_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
%This selects current working folder;
%Selected folder should contain all the data.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Current_Folder = handles.Current_Folder;
Selected_Folder = uigetdir(Current_Folder,'Select Folder with All Data');
Manual_selection = get(handles.checkbox1,'Value');
if (Selected_Folder ~= 0)
    Current_Folder = Selected_Folder;
    
    %analyze current folder to make sure it contains data
    Data_List = analyze_folder(Current_Folder, Manual_selection);
    handles.Data_List = Data_List;
    
    got_data = ~isempty(Data_List);
    if got_data
        selection = cell(1,length(Data_List));
        for i = 1 : length(Data_List)
            selection{i} = Data_List{i}.Description;
        end;
        set(handles.popupmenu1,'String',selection);
        set(handles.popupmenu1,'Value',1);
        set(handles.text2,'String',Current_Folder);
        set(handles.uipanel1,'Visible','on');
        
        
    else
        set(handles.text2,'String',['NO DATA Found in: ' Current_Folder]);
        set(handles.popupmenu1,'String',['No data Found']);
        set(handles.uipanel1,'Visible','off');
        
    end;
    
else
    if (strcmp(Current_Folder,''))
        set(handles.text2,'String','No Folder Selected! Please select folder...');
        set(handles.uipanel1,'Visible','off');
    end;
end;


handles.Current_Folder = Current_Folder;
guidata(hObject, handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Plot selection button
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles.All_Figures;
if ~isempty(h)
   for i = 1:length(h)
       try close(h{i}); catch exception;
       end;
   end;
end;
set(handles.listbox1,'String',{});
set(handles.listbox1,'Value',1);
selection = get(handles.popupmenu1,'Value');
Current_Data = handles.Data_List{selection};
[code, output, figures] = plot_all_data(Current_Data.PEG,Current_Data.KCL,...
    [handles.Current_Folder '\' Current_Data.Description ], Current_Data.Three_Bundle);
if (code ~= 0) 
    msgbox('Problem Plotting Data!' ,'Error!');
else
    set(handles.listbox1,'String',output);
    handles.All_Figures = figures;
    
end;
guidata(hObject,handles);



% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
