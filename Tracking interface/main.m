function varargout = main(varargin)
% MAIN M-file for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 21-Jan-2015 15:57:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_OpeningFcn, ...
    'gui_OutputFcn',  @main_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end


if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end;


% End initialization code - DO NOT EDIT


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
clc;


% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
%     Load file button
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Determine if this is the first time folder is opened;
existing_path = '';
try
    existing_path = handles.Properties.path;
catch any_exception;
end;

%now the program goes to the default path folder or to the existing path if
%image was already opened;

[FileName,PathName,FilterIndex] = uigetfile('.tif','Select TIFF file...',existing_path);


if (FileName ~= 0)
    
    clear_messages(handles);
    set_initial_gui_state(handles);
    
    msg = 'Loading file: ';
    update_messages(handles, [msg FileName ' ']);
    
    %Properties store run options,images, etc
    Properties.path = PathName;
    Properties.name = [PathName FileName];
    Properties.namec = '';
    Properties.id = '';
    Properties.day = '';
    
    Properties.threshold = 0.3;
    
    Properties.info = imfinfo(Properties.name);
    Properties.number_of_images = numel(Properties.info);
    Properties.current_image = 1;
    
    
    %if (matlabpool('size') == 0) - Due to removal of matlabpool in future
    %releases of Matlab
    
    
    if isempty(gcp('nocreate'))
        Properties.parallel_started = false;
    else
        Properties.parallel_started = true;
    end;
    
    
    %make panel actions with images visible;
    set(handles.uipanel1, 'Visible', 'On');
    
    %setup panel for manipulation with loaded images;
    set(handles.uipanel4, 'Title', Properties.name);
    set(handles.text5, 'String', ['Total Images: ' num2str(Properties.number_of_images)]);
    set(handles.text7, 'String', ['Current Image: ' num2str(Properties.current_image)]);
    %following part updates slider properties according to loaded images
    set(handles.slider4,'Min',1);
    set(handles.slider4,'Max', Properties.number_of_images);
    step_size = 1/(Properties.number_of_images - 1);
    set(handles.slider4, 'SliderStep', [step_size, step_size*10]);
    set(handles.slider4,'Value',1);
    set(handles.uipanel4, 'Visible', 'on');
    
    Properties.images = double(imread(Properties.name,Properties.current_image,'Info',Properties.info));
    
    image_handle = imagesc(Properties.images);
    colormap('gray');axis image; title (['Original Image # ' num2str(Properties.current_image)]);
    %plot(handles.axes1,image_handle);
    
    
    handles.Properties = Properties;
    guidata(hObject,handles);
    
    update_messages(handles, ['Total of ' num2str(Properties.number_of_images) ' images loaded.']);
    
end;


% --- Executes on button press in pushbutton2.
% --- Initialize parallel mode;
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_messages(handles, 'Initializing parallel mode...')

try
   poolobj = parpool;
catch exception
    update_messages(handles, [exception.message '\n']);
end;

%workers = matlabpool('size'); - due to removal of matlabpool in futute
%releases of Matlab;

workers = poolobj.NumWorkers;
update_messages(handles,...
    ['Ready for parallel computation with ' num2str(workers)  ' workers.']);
set(hObject,'Enable','off');

handles.Properties.parallel_started = true;
guidata(hObject,handles);



% --- Executes on button press in pushbutton3 - Track image
% --- Only tracks current image to see if the code works...
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%inverted and bypass-filtered images
image = handles.Properties.bp_images;
image_number = handles.Properties.current_image;

max_intensity = max(max(image));
threshold = handles.Properties.threshold * max_intensity;

%find peaks to pixel accuracy
pk = pkfnd(image,threshold,handles.Properties.max_bpass);
size_of_pk = size(pk);
num_of_peaks = size_of_pk(1);

if (num_of_peaks > 0)
    
    update_messages(handles,['Found ' num2str(num_of_peaks) ' particles.']);
    
    %get coordinates to subpixel resolution
    cnt=cntrd(image,pk,handles.Properties.max_bpass,0);
    
    
    %plot them the on the image
    set(handles.axes1,'NextPlot','add');
    for i = 1:num_of_peaks
        
        plot(pk(i,1),pk(i,2),'or');
        text(pk(i,1) + 10, pk(i,2), ['\leftarrow' num2str(i)],'Color','g');
        update_messages(handles, ['Particle # ' num2str(i) '...']);
        update_messages(handles,['X' num2str(i) ' = ' num2str(pk(i,1))...
            ' Y' num2str(i) ' = ' num2str(pk(i,2))]);
        plot(cnt(i,1),cnt(i,2),'*g');
        update_messages(handles,['Subpixel Results:  X' num2str(i) ' = ' num2str(cnt(i,1))...
            ' Y' num2str(i) ' = ' num2str(cnt(i,2))]);
        
    end;
    title(['Tracked Image # ' num2str(image_number)]);
else
    update_messages(handles,'No peaks found!');
end;

set(handles.axes1,'NextPlot','replace');


%make button7 - Track All enabled
set(handles.pushbutton7,'Enable','on');


% --- Executes on button press in pushbutton4
%     Bandpass filtering of the image
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = handles.Properties.images;
image_number = handles.Properties.current_image;

%min bandpass value
min_bpass = round(get(handles.slider1,'Value'));

%bandpass values should be odd
if (mod(min_bpass,2) == 0 )
    min_bpass = min_bpass + 1;
end;
handles.Properties.min_bpass = min_bpass;

%max bandpass value
max_bpass = round(get(handles.slider2,'Value'));

%bandpass values should be odd
if (mod(max_bpass,2) == 0 )
    max_bpass = max_bpass + 1;
end;


handles.Properties.max_bpass = max_bpass;

bp_images = bpass(image,min_bpass,max_bpass);

%new images are bandpass-filtered
handles.Properties.bp_images = bp_images;

guidata(hObject,handles);

image_handle = imagesc(bp_images);
colormap('gray');axis image; title (['Bandpass Filtered Image # ' num2str(image_number)]);

update_messages(handles,'Image has been filtered with bandpass');
update_messages(handles,['Min: ' num2str(min_bpass) '; Max: ' num2str(max_bpass) '.']);

set(handles.pushbutton3,'Enable','on');




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


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
% Hint: place code in OpeningFcn to populate axes1


% --- Executes when user attempts to close figure1.
% --- we close matlabpool workers here
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
 
 % matlabpool will be removed from Matlab future release; to comply, the
 % code  was modified
 %
 % try matlabpool close;
 % catch exception;
 % end;

try delete(gcp('nocreate'));
catch exception;
end;



function update_messages(handles, message)
%Updates messages in ListBox

current_msg = get(handles.listbox1, 'String');
current_msg{length(current_msg)+1} = message;
set(handles.listbox1, 'String', current_msg);
set(handles.listbox1, 'Value', length(current_msg));
drawnow;

function clear_messages(handles)
%Clears all messages in ListBox

set(handles.listbox1, 'String', '');
set(handles.listbox1, 'String', '');
set(handles.listbox1, 'Value', 1);

function set_initial_gui_state(handles)
%sets gui to initial state when new image is loaded
%except for Parallel Workers - this state remains the same when 
%the new image is loaded

set(handles.pushbutton3,'Enable','off');
set(handles.pushbutton4,'Enable','off');
set(handles.pushbutton7,'Enable','off');
set(handles.uipanel2,'Visible','off'); %panel for actions with tracked beads

set(handles.uipanel4,'Visible','off'); %panel with loaded images

set(handles.pushbutton13,'Enable','off');%Fourier Transform Separation button
set(handles.pushbutton14,'Enable','off');%Save Separation button
set(handles.pushbutton16,'Enable','off');%the Save Coordinate button

function scale = camera_scale(handles)
% gets scale in micrometers per pixel of the camera
%

phantom_cam = get(handles.radiobutton1,'Value');
phantom_zoomin = get(handles.radiobutton2,'Value');
phantom_zoomin_1point5x = get(handles.radiobutton8,'Value');
phantom_zoomin_1point5_and_2point5 = get(handles.radiobutton10,'Value');
andor_clara = get(handles.radiobutton5,'Value');
andor_ikon = get(handles.radiobutton6,'Value');
andor_ikon_zoomout = get(handles.radiobutton7,'Value');
eyeball = get(handles.radiobutton10,'Value');

%set magnification to 1(default, no scale);
scale = 1;

if phantom_cam
    scale = 0.115;
end;

if phantom_zoomin
    scale = 0.115/2.5;
end;

if phantom_zoomin_1point5x
    scale = 0.115/1.5;
end;

if phantom_zoomin_1point5_and_2point5
    scale = 0.115/(1.5*2.5);
end;


if andor_clara
    scale = 0.0645;
end;

if andor_ikon
    scale = 0.13;
end;

if andor_ikon_zoomout
    scale = 0.130/0.55;
end;

if eyeball
    msg_string = {'Ask Steve to measure distance between beads. That''s his idea...';...
        ''; 'Meanwhile, I will continue in pixels.'};
    h_msg = msgbox(msg_string,'Wanna Eyeball Separation?');
    uiwait(h_msg);
    set(handles.radiobutton10,'Value',0);
    set(handles.radiobutton4,'Value',1);
end;



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
msg_string = ['Min bpass value: ', num2str(round(get(hObject,'Value')))];
%update_messages(handles, msg_string);
set(handles.text3,'String',msg_string);



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider2 movement.
% --- Setting max bandpass value
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%get the reading
msg_string = ['Max bpass value: ', num2str(round(get(hObject,'Value')))];
%update_messages(handles, msg_string);
set(handles.text4,'String',msg_string);



% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton6
%     Invert image
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.Properties.images;
image_number = handles.Properties.current_image;
bit_depth = handles.Properties.info.BitDepth;

inverted = invert_image(image,bit_depth);
handles.Properties.images = inverted;
handles.Properties.original = image;
guidata(hObject,handles);

image_handle = imagesc(inverted);
colormap('gray');axis image; title (['Inverted Image # ' num2str(image_number)]);

update_messages(handles,'Image has been inverted');
set(handles.pushbutton4,'Enable','on');


% --- Executes during object creation, after setting all properties.
function text3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton7
% --- Track All Images button...
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if (~handles.Properties.parallel_started)
    
    update_messages(handles,'Parallel Mode Required...');
    
else
    
    set(handles.uipanel1,'Visible','off'); drawnow expose;
    
    update_messages(handles,'Please Wait, Tracking Beads ...');
    
    
    h = handles.Properties;
    bit_depth = h.info.BitDepth;
    
    %calculate time to track
    start_time = tic;
    
    handles.Properties.tracked_coords = beadtrack_parallel(h.name,...
        bit_depth, h.min_bpass, h.max_bpass);
    
    
    elapsed_time = toc(start_time);
    
    %number of tracked particles
    num_tracked = max(handles.Properties.tracked_coords(:,4));
    handles.Properties.num_tracked = num_tracked;
    
    update_messages(handles,['Tracked ' num2str(num_tracked) ' beads in ' num2str(elapsed_time) ' s.'] );
    
    
    %make visible panel with additional buttons
    if num_tracked > 0
        
        update_messages(handles,['Time per frame: ', num2str(elapsed_time/h.number_of_images) ' s.']);
        %populate drop-down menu for individual bead's coordinates
        %populate drop-down menus for separation
        bead_list = cell(num_tracked,1);
        bead_coord_list_x = cell(num_tracked,1);
        bead_coord_list_y = cell(num_tracked,1);
        
        for i = 1:num_tracked
            bead_list{i} = ['Bead ' num2str(i)];
            bead_coord_list_x{i} = ['X' num2str(i)];
            bead_coord_list_y{i} = ['Y' num2str(i)];
        end;
        bead_coord_list = [bead_coord_list_x; bead_coord_list_y];
        
        set(handles.popupmenu1,'String', bead_list);
        set(handles.popupmenu2,'String', bead_list);
        set(handles.popupmenu3,'String', bead_list);
        set(handles.popupmenu5,'String', bead_coord_list);
        
        set(handles.popupmenu1,'Value', 1);
        set(handles.popupmenu2,'Value', 1);
        set(handles.popupmenu3,'Value', 1);
        set(handles.popupmenu5,'Value', 1);
        
        set(handles.uipanel2,'Visible','on');
    else
        set(handles.uipanel2,'Visible','off');
        update_messages(handles,'No beads detected, check filtering and try again!');
    end;
    
    set(handles.uipanel1,'Visible','on');
    guidata(hObject,handles);
    
end;


% --- Executes on button press in pushbutton8.
% --- Plot Position histogram for a given coordinate
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get what coordinate to plot...
value = get(handles.popupmenu5,'Value');
strings = get(handles.popupmenu5,'String');
selection = strings{value};
coordinate = selection(1);
bead_index = str2num(selection(2:end));

%coordinates are in all_coords
all_coords = handles.Properties.tracked_coords;

%select appropriate bead coords
bead_coords_indeces = find(all_coords(:,4) == bead_index);
bead_coords = all_coords(bead_coords_indeces,:);

%choose right coordinate
switch coordinate
    case'X'
        coord_index = 1;
    case 'Y'
        coord_index = 2;
end;

%determine scale
mag = camera_scale(handles);
if (mag == 1)
    units_string = 'pixels';
else
    units_string = ' \mum';
end;

%this is selected coordinate scaled appropriately
coord = mag * bead_coords(:,coord_index);

%plot histogram
figure(3);
%find out optimal number of bins using appropriate formula
[n_bins, bin_size] = get_optimal_binning(coord);
update_messages(handles,['Plot histograms with ' num2str(n_bins) ' bins.']);
hist(coord,n_bins);

%evaluate standard deviation of bead (in nanometers). Coord has coordinates
%in microns (all default length variables are in microns), so multiply by 
%1000 to get nm.
bead_1_x_std = std(coord)*1000;


%evaluate trap stiffness in given direction from standard deviation
% p(x) ~ exp(-kx^2/2(kBT)) ~ exp(-x^2/2s^2)
%therefore k = kBT/s^2
%use kBT = 4.14 pN * nm at room T
%stiffness is in pN/nm
stiffness = 4.14/(bead_1_x_std)^2;

h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','b');

%create title for histogram
title_string = ...
    {['\fontsize{12}' selection ' histogram'],[], ['\fontsize{10} \delta(' selection(1) '_' num2str(bead_index) ') = '...
    num2str(bead_1_x_std,3) ' nm';], ['Calcultate trap stiffness k = '...
    num2str(stiffness,3) ' pN/nm']};
title(title_string);
xlabel(['Bead location, ' units_string]);
%plot time-series;
figure(4);
plot(coord,'.r');

%save coordinate vs frame plot in GUI handles for possible save
handles.Properties.selected_coordinate = coord;
guidata(hObject,handles);
%enable the button to save coordinate time-series;
set(handles.pushbutton16,'Enable','on');





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


% --- Executes on button press in pushbutton9.
% --- Plot separation time-series
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%disable button
%set(hObject,'Enable','off');

%determine scale
mag = camera_scale(handles);

if (mag == 1)
    units_string = 'pixels';
else
    units_string = ' \mum';
end;

%get values of the two pop-up menus
bead_1 = get(handles.popupmenu2,'Value');
bead_2 = get(handles.popupmenu3,'Value');

%if same bead is selected, plot bead distance from (0,0)
if (eq(bead_1,bead_2))
    update_messages(handles,'Can not plot separation! Same bead is selected.');
    update_messages(handles,'Plotting bead time series instead...');
    
    %get coordinates for the first
    all_coords = handles.Properties.tracked_coords;
    bead_coords_indeces = find(all_coords(:,4) == bead_1);
    bead_coords_1 = all_coords(bead_coords_indeces,:);
    bead_coords_1(:,4) = 1;
    
    %second bead is now just a reference point at 0
    bead_coords_2 = bead_coords_1;
    for i = 1:length(bead_coords_2)
        bead_coords_2(i,[1,2]) = [0, 0];
        bead_coords_2(i,4) = 2;
    end;
    
    %generate title of the plot
    title_msg = ['Bead ', num2str(bead_1),' trajectory.'];
    
else
    update_messages(handles,'Plotting separation time series.');
    %get coordinates for the first and second bead
    all_coords = handles.Properties.tracked_coords;
    bead_coords_indeces = find(all_coords(:,4) == bead_1);
    bead_coords_1 = all_coords(bead_coords_indeces,:);
    
    %since we are not limited by two beads, rename last index to 1 for
    %the first bead
    bead_coords_1(:,4) = 1;
    
    
    bead_coords_indeces = find(all_coords(:,4) == bead_2);
    bead_coords_2 = all_coords(bead_coords_indeces,:);
    
    %index must be two here in order for this to work
    bead_coords_2(:,4) = 2;
    
    %generate title of the plot
    title_msg = ['Separation between beads ', num2str(bead_1),' and ', num2str(bead_2),'.'];
end;

%make one vector from bead_coords_1,2
bead_coords = [bead_coords_1; bead_coords_2];

%use SortBeads to get trajectories
trajectories = SortBeads(bead_coords);

%calculate separation
separation = magR(trajectories,1,3) * mag;

figure(4);
p = plot(separation);
axis tight;
set(p,'LineStyle','none','Marker','.','Color','red');
title(title_msg);
xlabel('Frame #');
ylabel(['Distance,' units_string]);

%save given trajectory in handles
handles.Properties.separation = separation;
guidata(hObject,handles);

%show FT button...
set(handles.pushbutton13,'Enable','on');

%show save separation button
set(handles.pushbutton14,'Enable','on');




% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
% --- Check for Randomness
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bead_index = get(handles.popupmenu1,'Value');

all_coords = handles.Properties.tracked_coords;
bead_coords_indeces = find(all_coords(:,4) == bead_index);
bead_coords = all_coords(bead_coords_indeces,:);

f1 = figure(1);
hist(mod(bead_coords(:,1), 1));
title(['X - coordinate remainders for bead ' num2str(bead_index)]);
x_rems = hist(mod(bead_coords(:,1), 1));

f2 = figure(2);
y_rems = hist(mod(bead_coords(:,2), 1));
hist(mod(bead_coords(:,2), 1));
title(['Y - coordinate remainders for bead ' num2str(bead_index)]);

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton13.
% --- Fourier Transform trajectory from time domain
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_messages(handles,'Plotting FT of trajectory...');

%first, get the trajecotry
separation = handles.Properties.separation;

update_messages(handles,'Attempting FT...');

[fc, D] = fourier_transform_routine(separation);

update_messages(handles,['Corner frequency estimated: ' num2str(fc) ' Hz']);
update_messages(handles,['Diffusion coeff est: ' num2str(D) 'mum/s^2']);




% --- Executes on button press in pushbutton14.
% --- Save separation
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get relevant separation that needs to be saved
separation = handles.Properties.separation;

%get directory name where images are stored
%and image filename;
filename = handles.Properties.name;
last_char_index = find(filename == '\', 1, 'last' );
directory = filename(1:last_char_index);
name = filename(last_char_index + 1 : end);


%open save dialog window;
save_dlg(separation, directory, name);


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try filename = handles.Properties.name;
catch exception;
    filename = userpath;
end;

last_char_index = find(filename == '\', 1, 'last' );
directory = filename(1:last_char_index);

load_dlg(0,directory);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Properties = handles.Properties;
Properties.current_image = round(get(hObject,'Value'));
set(hObject,'Value',Properties.current_image);
set(handles.text7, 'String', ['Current Image: ' num2str(Properties.current_image)]);

%plot current image
Properties.images = double(imread(Properties.name,Properties.current_image,'Info',Properties.info));
image_handle = imagesc(Properties.images);
colormap('gray');axis image; title (['Original Image # ' num2str(Properties.current_image)]);

handles.Properties = Properties;
guidata(hObject,handles);





% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton16.
%Saves the selected coordinate separation in .mat file
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coord = handles.Properties.selected_coordinate;
filename = handles.Properties.name;
last_char_index = find(filename == '\', 1, 'last' );
directory = filename(1:last_char_index);
name = filename(last_char_index + 1 : end);

save_dlg(coord, directory, name);
