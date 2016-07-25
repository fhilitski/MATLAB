function varargout = FilamentTrackerII(varargin)
% FILAMENTTRACKERII M-file for FilamentTrackerII.fig
%      FILAMENTTRACKERII, by itself, creates a new FILAMENTTRACKERII or raises the existing
%      singleton*.
%
%      H = FILAMENTTRACKERII returns the handle to a new FILAMENTTRACKERII or the handle to
%      the existing singleton*.
%
%      FILAMENTTRACKERII('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILAMENTTRACKERII.M with the given input arguments.
%
%      FILAMENTTRACKERII('Property','Value',...) creates a new FILAMENTTRACKERII or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FilamentTrackerII_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FilamentTrackerII_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FilamentTrackerII

% Last Modified by GUIDE v2.5 19-Mar-2015 18:38:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FilamentTrackerII_OpeningFcn, ...
    'gui_OutputFcn',  @FilamentTrackerII_OutputFcn, ...
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

% --- Executes just before FilamentTrackerII is made visible.
function FilamentTrackerII_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FilamentTrackerII (see VARARGIN)

% Choose default command line output for FilamentTrackerII
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FilamentTrackerII wait for user response (see UIRESUME)
% uiwait(handles.figure1);
clc;
%Load Image File
tmp_path = 'D:\Data - MT Sliding and Friction\2015\*.*';

name = 0;
frame_index = 1;
handles.frame_index = frame_index;

%ask user for the image file first
[name,file_path] = uigetfile(tmp_path,'Select Image file');
if name == 0
    close(handles.figure1);
else
    handles.Imgname = name;
    handles.Imgpathname = file_path;
        
    handles.ImInfo=imfinfo([handles.Imgpathname,handles.Imgname]);
    handles.Nmax=length(handles.ImInfo);
    set(handles.slider4,'Max',handles.Nmax);
    slider_step = 1/handles.Nmax;
    set(handles.slider4,'SliderStep',[slider_step, slider_step]);
    set(handles.slider4,'Value',frame_index);
    set(handles.edit5,'String',num2str(frame_index));
    
    handles.Img=imread([handles.Imgpathname,'\',handles.Imgname],frame_index);
    
    %Load Coordinate File
    %assume it is in the same folder as images
    coord_file_path = handles.Imgpathname;
    %name is now something like 'image.ext', for example 'images.tif'
    %we need to extract the actual filename, and drop extension.
    [tmp_path, file_name, last_extension] = fileparts(name);
    %now file_name is just 'images';
    %traditionally, we store tracked skeletons in '.skl.txt' file
    name = [file_name '.skl.txt'];
    try 
        handles.coords = load([file_path,'',name]);
    catch exception_on_load
        disp(['Could not find ' name]);
        name = 0;
        while (name == 0)
            [name,file_path] = ...
                uigetfile([coord_file_path '\*.*'],'Select Coordinate File');
        end;
        handles.coords = load([file_path,'\',name]);
    end;
     
    msg_string = ['Loaded: ' file_path name];
    set(handles.text5,'String',msg_string);
    guidata(hObject, handles);
    %Find Filament Coordinates
    handles = get_contour(hObject, eventdata, handles);
    
    % Get relevant parameters: smoothing, width, threshold
    % for the first time.
    handles.s = get(handles.slider1,'Value');
    handles.w = get(handles.slider2,'Value');
    handles.minpeak = get(handles.slider3,'Value');
    handles.imgcount=1;
    
    
    handles = get_intensity_profiles(hObject, eventdata, handles);
    plot_all_images(hObject, eventdata, handles,frame_index);
    update_profile_plots(hObject,eventdata,handles);
end;
set_initial_gui_state(hObject, eventdata, handles);

%obtain intensity profile for image in handles structure;
function h = get_intensity_profiles(hObject, eventdata, handles)

handles.IntPr = IntProfile_ANDY(handles.contour,handles.Img,handles.w);
handles.IntPrS = smooth(handles.IntPr,handles.s);
handles.dIntPrS=handles.IntPrS(2:end)-handles.IntPrS(1:end-1);

guidata(hObject,handles);
h = handles;



%obtains contour for image and stores it in handles
function h =  get_contour(hObject, eventdata, handles)
%get sekeleton
handles.f = find(handles.coords(:,1) == handles.frame_index);
handles.f1 = handles.f(1);
handles.L = length(handles.f); %Length of the filament
handles.contour = handles.coords(handles.f,2:3);
%handles.contour(:,1)=handles.contour(:,1)+1;
[handles.CL,handles.contour]=ContourLength(handles.contour);

%save skeleton
guidata(hObject,handles);
h = handles;


%plots selected images and skeleton
function plot_all_images(hObject, eventdata, handles, image_index)

handles.Img=imread([handles.Imgpathname,handles.Imgname],image_index);

axes(handles.Image);
imshow(handles.Img,[min(min(handles.Img)) max(max(handles.Img))]);
%find appropriate skeleton
handles.f = find(handles.coords(:,1) == image_index);
handles.f1 = handles.f(1);
handles.L = length(handles.f); %Length of the filament
handles.contour = handles.coords(handles.f,2:3);
%handles.contour(:,1)=handles.contour(:,1)+1;

[handles.CL,handles.contour]=ContourLength(handles.contour);
hold on;
plot(handles.contour(:,1),handles.contour(:,2),'LineWidth',2)



%Plots intensity profile, its derivative and threshold with set parameters.
function update_profile_plots(hObject, eventdata, handles)

%calculate Intensity Threshold
minInt=(1-handles.minpeak)*max(handles.IntPrS);
maxInt=(1+handles.minpeak)*max(handles.IntPrS);

%plor smoothed profile
axes(handles.Profile);
plot(handles.CL,handles.IntPrS);
hold on;
%plot smoothed derivative
plot(handles.CL(1:end-1),abs(handles.dIntPrS),'r');
%plot threshold lines
handles.h_top = plot(handles.CL,maxInt,'-b');
handles.h_bottom = plot(handles.CL,minInt,'-b');
%adjust axes limits
axis tight;
y = ylim;
y(2) = y(2) + 10;
ylim(y);
hold off;






% --- Outputs from this function are returned to the command line.
function varargout = FilamentTrackerII_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;

% --- Executes on slider1 movement
% --- slider regulates Intensity Smoothing
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = round(get(hObject, 'Value'));
set(handles.Smooth, 'String', num2str(s));
handles.s = s;
set(hObject,'Value',s);

%recalculate smoothed profiles
handles.IntPrS = smooth(handles.IntPr,handles.s);
handles.dIntPrS=handles.IntPrS(2:end)-handles.IntPrS(1:end-1);
%save smoothed profiles
guidata(hObject,handles);
update_profile_plots(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function Smooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement
% --- Changing "intensity width"
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
w = round(get(hObject,'Value'));
set(handles.Width, 'String', num2str(w));
set(hObject,'Value',w);
handles.w = w;

%update all profiles
handles.IntPr = IntProfile_ANDY(handles.contour,handles.Img,w);
handles.IntPrS = smooth(handles.IntPr,handles.s);
handles.dIntPrS = handles.IntPrS(2:end)-handles.IntPrS(1:end-1);
%save updated profiles and plot them
guidata(hObject, handles);
update_profile_plots(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function Width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function PeakH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeakH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ???
% There is no such object at this point
% ???
function MinPeakW_Callback(hObject, eventdata, handles)
% hObject    handle to MinPeakW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinPeakW as text
%        str2double(get(hObject,'String')) returns contents of MinPeakW as a double
guidata(hObject, handles);
handles.minPwidth=str2double(get(handles.MinPeakW,'String'));
guidata(hObject, handles);

% ???
% There is no such object at this point
% ???
% --- Executes during object creation, after setting all properties.
function MinPeakW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinPeakW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Find Filament.
% ---
function FindPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to FindPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
threshold = handles.minpeak;
minInt=(1-threshold)*max(handles.IntPrS);
maxInt=(1+threshold)*max(handles.IntPrS);

%get overlap contour
f=find(handles.IntPrS<=maxInt & handles.IntPrS>=minInt);
CLf = handles.CL(f);
IntPf=handles.IntPrS(f);

%midpoint of the bundle
midp=(CLf(end)+CLf(1))/2;

%average filament intensity
minI=mean(IntPf);

%plot filament on intensity profile plot
axes(handles.Profile);
%h_overlap_plot = figure;
hold on;
%plot(h_overlap_plot,handles.CL,handles.IntPrS);
plot(CLf,IntPf,'go');plot(midp,minI,'ro','MarkerFaceColor','red');
hold off;

set(handles.Track,'Enable','on');
guidata(hObject, handles);


function[peakdata]=GetPeaks(CL,IntP,minP,minW)
%function to find peaks in intensity profile
[peakdata]=mspeaks(CL',IntP,'HeightFilter',minP,'FWHHFilter', minW);


% --- Executes on button press in Track.
function Track_Callback(hObject, eventdata, handles)
% hObject    handle to Track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = handles.Imgpathname;
imgname = handles.Imgname;
coords = handles.coords;
width = handles.w;
smoothing = handles.s;

%total = handles.Nmax;
total = max(coords(:,1));


x=zeros(total, 4);
%h_tmp_fig = figure;
count = handles.imgcount;

for count = 1:total
    
    %Read Image
    Img = imread([path,'\',imgname],count);
    
    %Find Filament Coordinates
    f = find(coords(:,1) == count);
    
    %in case the filament is not found in a given frame.
    if isempty(f)
        disp(['Filament not found, frame #' num2str(count)]);
        x(count,:) = [count,0,0,0];
        continue
    end;
    
    %also, there are problems with really short filaments
    %define a threshold amount of data-points necessary to track a single
    %filament in a frame
    points_threshold = 25;
    if (length(f) < points_threshold)
        disp(['Filament is too short for tracking in frame #', num2str(count)]);
        x(count,:) = [count,0,0,0];
        continue
    end;
    f1 = f(1);
    
    %Get Contour
    contour = coords(f,2:3);
    %handles.contour(:,1) = handles.contour(:,1) + 1;
    
    [CL, contour] = ContourLength(contour);
    
    %     if CL(end) > 1.2*filament_length
    %         count = count+1;
    %         continue
    %     end
    %
    %Check that the order is correct
    %     if count > 1
    %        % dS1=sqrt((contour(1,1) - contourp(1,1))^2+(contour(1,2)-contourp(1,2))^2);
    %        % dS2=sqrt((contour(1,1) - contourp(end,1))^2+(contour(1,2)-contourp(end,2))^2);
    %         % #ok<COLND>
    %         if dS2<dS1
    %             contour=contour(end:-1:1,:);
    %             [CL,contour]=ContourLength(contour);
    %         end
    %     end
    
    
    contourp = contour;
    
    %Get Intensity Profile
    IntPr=IntProfile_ANDY(contour, Img, width);
    
    %Smooth Data and take derivative
    IntPrS=smooth(IntPr,smoothing);
    %handles.dIntPrS=handles.IntPrS(2:end)-handles.IntPrS(1:end-1);
    
    
    %Find Filament
    minInt=(1-handles.minpeak)*max(IntPrS);
    maxInt=(1+handles.minpeak)*max(IntPrS);
    f=find( IntPrS <= maxInt & IntPrS>=minInt);
    CLf = CL(f);
    IntPf = IntPrS(f);
    midp=(CLf(end)+CLf(1))/2;
    minI = mean(IntPf);
    
    x(count,:) = [count,midp,CL(end),CLf(end)-CLf(1)];
    %
    %     figure(h_tmp_fig);
    %     hold off;
    %     subplot(2,1,1);
    %     imshow(Img,[min(min(Img)) max(max(Img))]);
    %     hold on;
    %     plot(contour(:,1),contour(:,2));
    %     subplot(2,1,2);
    %     hold off;
    %     plot(CL,IntPrS);
    %     hold on;
    %     plot(CLf,IntPf,'go');
    %     plot(midp,minI,'ro','MarkerFaceColor','red');
    %     title(['Image: ' num2str(count) ' out of ' num2str(total)]);
    if (mod(count, 1) == 0)
        disp(['Image: ' num2str(count) ' out of ' num2str(total)]);
    end;
end;
disp('Done!');
msgbox('Tracking is complete!','FilamentTracker');
handles.x = x;

set(handles.ExportData,'Enable','on');
set(handles.pushbutton5,'Enable','on');
set(handles.pushbutton6,'Enable','on');
guidata(hObject, handles);


% --- Executes on button press in Export Diffusion Data.
function ExportData_Callback(hObject, eventdata, handles)
% hObject    handle to ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%this script analyzes and exports bundle diffusion data, position of the
%bundle and plots MSDs. It is important that the total length of the
%filament remains constant (i.e. there is no elongation or shrinkage of
%microtubules. Datapoints where total measured length is different
%significantly are thrown out by process_filament_data script, as the
%filament is likely moved out of focus.
%To track growth of the filaments, use Export Length Data button.

%handles.x contains  [count,midp,CL(end)]
assignin('base', 'Xt', handles.x);
x = handles.x;
xvst = x(:,2);
L = x(:,3);
assignin('base', 'x', xvst);
assignin('base', 'L', L);
axes(handles.Trajectory);
hold off;
plot(x(:,1),x(:,2),'.r');
hold on;
plot(x(:,1),x(:,3),'.b');
axis tight;
title('Raw Tracking Results');

[time, msd, error, new_x, new_t] = process_filament_data(handles.x,0.5,0.1,1);

fname = [handles.Imgname '.mat'];
fname = [handles.Imgpathname fname];
save(fname, 'new_x','new_t','x','time', 'msd', 'error');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider3 movement.
% --- This slider adjusts intensity threshold.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%Get intensity threshold value
peakH = get(hObject,'Value');
set(handles.PeakH, 'String', num2str(peakH));
handles.minpeak = peakH;

%save new threshold value
guidata(hObject, handles);
update_profile_plots(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
selected_frame = round(get(hObject,'Value'));
set(hObject,'Value',selected_frame);
set(handles.edit5,'String',num2str(selected_frame));
handles.frame_index = selected_frame;

%read selected image
handles.Img=imread([handles.Imgpathname,'\',handles.Imgname],selected_frame);
guidata(hObject,handles);
%Find Filament Coordinates
handles = get_contour(hObject, eventdata, handles);
handles = get_intensity_profiles(hObject, eventdata, handles);
plot_all_images(hObject, eventdata, handles,selected_frame);
update_profile_plots(hObject,eventdata,handles);
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


% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on export length data press(Export Length Data)
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Exports length of the full filament
% Threshold and width values are not relevant

%handles.x contains  [count,midp,CL(end)]
%assignin('base', 'Xt', handles.x);
x = handles.x;
t = x(:,1); %frame number
xvst = x(:,2); %midpoint position; not relevant for length measurement
L = x(:,3); %length of the filament
L_segment = x(:,4); %length of the thresholded segment


axes(handles.Trajectory);
hold off;
plot(t,L,'.r');
axis tight;
title('Raw Filament Length');
xlabel('frames');
ylabel('pixels');

%ask for the time and length-scale
user_input = {};
while isempty(user_input)
    prompt = {'Pixel size, microns'; 'Frame time, ms'};
    def_values = {'0.13';'0.118'};
    user_input = inputdlg(prompt,'Set the scale',1,def_values);
end;
pixel_size = str2num(user_input{1}); %microns per pixel
time_scale = str2num(user_input{2}); %ms per frame

[time_frames, l_total, time_s, length_um, segment_um] = process_length_data(handles.x,time_scale,pixel_size);

%create filename for the .mat file to save the length data
fname = [handles.Imgname '.length.mat'];
fname = [handles.Imgpathname fname];
%save relevant variables into the mat file
save(fname, 'time_frames','l_total','time_s','length_um','segment_um');
guidata(hObject, handles);
cftool(time_s, length_um);


% --- Executes on button press in pushbutton6 (Export Overlay Images).
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = handles.Imgpathname;
imgname = handles.Imgname;
coords = handles.coords;
image_index = 1;
img=imread([pathname,imgname],image_index);

%find appropriate skeleton
f = find(coords(:,1) == image_index);
f1 = handles.f(1);
L = length(f); %Length of the filament
contour = coords(f,2:3);

%plot image and contour
% h_img_fig = figure;
% h_subplot_left = subplot(1,2,1);
% imshow(img,[min(min(img)) max(max(img))]);
% h_subplot_right = subplot(1,2,2);
% imshow(img,[min(min(img)) max(max(img))]);
% hold on;
% plot(contour(:,1),contour(:,2),'r','LineWidth',2);


%create an overlay image
overlay = addContour(img,contour,[0,0],false);
h_overlay_fig = figure;
imshow(overlay);

%determine its size
[m,n,l] = size(overlay);
%pre-make frame array
total_frames = handles.Nmax;
all_frames = zeros(m,n,l,total_frames,'uint8');




for i = 1:total_frames
    disp(['Exporting frame: ', num2str(i)]);
    img=imread([pathname,imgname],i);
    
    f = find(coords(:,1) == i);
    f1 = handles.f(1);
    contour = coords(f,2:3);
    
    all_frames(:,:,:,i) = addContour(img,contour,[1,1],false);
    
end;
disp('Done exporting!');
video_writer = VideoWriter([pathname imgname],'Uncompressed AVI');
open(video_writer)
writeVideo(video_writer,all_frames)
close(video_writer);
msgbox('Video saved!','FilamentTracker');


%
function set_initial_gui_state(hObject, eventdata, handles)
%Function sets gui state after loadig the image.
set(handles.FindPeaks,'Enable','on');
set(handles.Track,'Enable','off');
set(handles.ExportData,'Enable','off');
set(handles.pushbutton5,'Enable','off');
set(handles.pushbutton6,'Enable','off');
guidata(hObject, handles);
