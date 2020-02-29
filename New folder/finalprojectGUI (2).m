function varargout = finalprojectGUI(varargin)
% FINALPROJECTGUI MATLAB code for finalprojectGUI.fig
%      FINALPROJECTGUI, by itself, creates a new FINALPROJECTGUI or raises the existing
%      singleton*.
%
%      H = FINALPROJECTGUI returns the handle to a new FINALPROJECTGUI or the handle to
%      the existing singleton*.
%
%      FINALPROJECTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINALPROJECTGUI.M with the given input arguments.
%
%      FINALPROJECTGUI('Property','Value',...) creates a new FINALPROJECTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before finalprojectGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to finalprojectGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help finalprojectGUI

% Last Modified by GUIDE v2.5 30-Mar-2016 05:15:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @finalprojectGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @finalprojectGUI_OutputFcn, ...
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


% --- Executes just before finalprojectGUI is made visible.
function finalprojectGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to finalprojectGUI (see VARARGIN)

% Choose default command line output for finalprojectGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes finalprojectGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = finalprojectGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = waitbar(0,'Loading Logo Database...');
logodict = 4;
mydata = cell(1, logodict);

for k = 1:logodict
  myfilename = sprintf('logo%d.jpg', k);
  eval(sprintf('I%d = imread(myfilename);',k));
  eval(sprintf('handles.logo%d= I%d',k,k))
end
handles.logodict = logodict;
close(h);
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hh = waitbar(0,'Please Wait...');
logodict = handles.logodict;
for i = 1:logodict
eval(sprintf('I%d= handles.logo%d;',i,i));
eval(sprintf('[a1 a2]= size(I%d);',i));
d(i)= a1*a2;
end

[m n] = min(d);
eval(sprintf('[a1 a2]= size(I%d);',n));

for i=1:logodict
eval(sprintf('I%d = imresize(I%d, [50 50]);',i,i));
figure('Visible','off');
eval(sprintf('h=imshow(I%d);',i));
saveas(h,sprintf('templogo_%d.jpg',i));
end


axes(handles.axes2);
fileFolder =  pwd;
dirOutput = dir(fullfile(fileFolder,'templogo_*.jpg'));
fileNames = {dirOutput.name}';

% for i=1:9
%     eval(sprintf('dirOutput.I%d = I%d;',i,i));
%     eval(sprintf('fileNames = {fileNames dirOutput.I%d};',i));
% end
montage(fileNames,[a1*9 a2*9]);
close(hh);




% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = waitbar(0,'Please wait...');
I = handles.Pic;




close(h)




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



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f1 = 1;
handles.f = f1; 
I= handles.Pic;
I = rgb2gray(I);
h = waitbar(0,'Please wait...');
[des1 loc1] = sift(I);
loc1 = loc1(:,1:2);
colormap('gray');
axes(handles.axes2);
imagesc(I);
hold on;
plot(loc1(:,2),loc1(:,1),'g+');
hold off;
handles.des1 = des1;
handles.loc1 = loc1;
    
close(h) 
% Choose default command line output for Example_GUI
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f2 = 2;
handles.f = f2; 
I = handles.Pic;
 I = rgb2gray(I);

h = waitbar(0,'Please wait...');

points = detectSURFFeatures(I);
    [features, valid_points] = extractFeatures(I, points);
    axes(handles.axes2);
imshow(I,'parent',handles.axes2); hold on;
    plot(valid_points.selectStrongest(size(features,1)),'showOrientation',true); hold off;
    des11 = features;
    loc11 = valid_points.Location;
    handles.des11 = des11; 
    handles.loc11 = loc11;

close(h) 
% Choose default command line output for Example_GUI
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f3 = 3;
handles.f = f3; 
I = handles.Pic;
h = waitbar(0,'Please wait...');

[x1 y1 v1]= harris(I);
loc111 = [x1 y1];
r1= ones(length(x1),1);
r1 = r1*10;
sif1 = find_sift(I,[x1 y1 r1],2);
des111 = sif1;
axes(handles.axes2);
colormap('gray');
imagesc(I);
hold on;
plot(x1,y1,'g+'); hold off
handles.des111 = des111;
handles.loc111 = loc111;  
close(h) 
% Choose default command line output for Example_GUI
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f4 = 4;
handles.f = f4; 
I = handles.Pic;
h = waitbar(0,'Please wait...');
corners   = detectFASTFeatures(rgb2gray(I));
     strongest = selectStrongest(corners, length(corners));
     [hog1, validPoints1, ptVis1] = extractHOGFeatures(I, strongest);
     axes(handles.axes2);
     imshow(I); hold on;
     plot(ptVis1, 'Color','green');
     
     loc1111 = validPoints1.Location;
     des1111 = hog1;
     handles.des1111 = des1111;
handles.loc1111 = loc1111; 
     
% title(sprintf('Location of %dx%d=%d SIFT patches of size = %dx%d' , options.deltay,options.deltax,options.deltay*options.deltax,options.patchsize,options.patchsize) ,'fontname' , 'times' , 'fontsize' , 13, 'fontweight','bold')
close(h)
% Choose default command line output for Example_GUI
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);



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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
logodict = handles.logodict;
h = waitbar(0,'Please wait...');
distRatio =get(handles.popupmenu4,'Value');
for i= 1:logodict
eval(sprintf('I%d = handles.logo%d;',i,i)); 
end
switch distRatio
    case 1
        rat=1;
    case 2
        rat=0.95;
    case 3
        rat=0.9;
    case 4
        rat=0.875;
    case 5
        rat=0.85;
    case 6
        rat=825;
    case 7
        rat=0.8;
    case 8
        rat=0.7;
    case 9
        rat=0.6;
    case 10
        rat=0.5;
end

f = handles.f;
% a = exist( 'handles.logo1', 'var');
%      if a== 0 
%         h = errordlg('Load Database');
%      end

I = handles.Pic;
kkk = handles.axes2;
if f ==1
des1 = handles.des1;
loc1 = handles.loc1;
loc1 = loc1(:,1:2);
pp = 'g+';
epsilon=30;  %dbscan search area
MinPts=5;     %dbscan minimun num of points for a cluster
axes(handles.axes2);
imagesc(I); hold on;
for i = 1:logodict
eval(sprintf('[des_logo%d, loc_logo%d] = sift(I%d);',i,i,i));
eval(sprintf('[matchLocx%d matchLocy%d num] = siftMatch(I, des1, loc1, des_logo%d, loc_logo%d, kkk , rat);',i,i,i,i));
eval(sprintf('handles.matchLocx%d = matchLocx%d;',i,i));
eval(sprintf('handles.matchLocy%d = matchLocy%d;',i,i));
eval(sprintf('idx%d=DBSCAN(matchLocx%d,epsilon,MinPts);',i,i));
eval(sprintf('handles.idx%d = idx%d;',i,i));
eval(sprintf('Class%d = unique(idx%d)',i,i));
eval(sprintf('Len%d = length(Class%d)',i,i));
eval(sprintf('handles.Len%d = Len%d;',i,i));
eval(sprintf('plot(matchLocx%d, matchLocy%d, pp);',i,i)); hold on;
end

elseif f ==2;
    des11 = handles.des11;
loc11 = handles.loc11;
I1 = rgb2gray(I1);
points = detectSURFFeatures(I1);
    [features, valid_points] = extractFeatures(I1, points);
    des2 = features;
    loc2 = valid_points.Location;   
    [matchLoc1 matchLoc2 num] = siftMatch(I, des11, loc11, des2, loc2, kkk , rat);
elseif f == 3;
    des111 = handles.des111;
loc111 = handles.loc111;
[x2 y2 v2]= harris(I1);
loc2 = [x2 y2];
r2= ones(length(x2),1);
r2 = r2*10;
sif2 = find_sift(I,[x2 y2 r2],2);
des2 = sif2;    
    [matchLoc1 matchLoc2 num] = siftMatch(I, des111, loc111, des2, loc2, kkk , rat);
elseif f == 4;
    des1111 = handles.des1111;
loc1111 = handles.loc1111;
corners   = detectFASTFeatures(rgb2gray(I1));
     strongest = selectStrongest(corners, length(corners));
     [hog2, validPoints2, ptVis2] = extractHOGFeatures(I1, strongest);
     loc2 = validPoints2.Location;
     des2 = hog2;
[matchLoc1 matchLoc2 num] = siftMatch(I, des1111, loc1111, des2, loc2, kkk, rat);
end

close(h);
% Choose default command line output for Example_GUI
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a= imaqhwinfo('winvideo',1); 
vid = videoinput('winvideo', 1, a.DefaultFormat);
vid.FramesPerTrigger = 1;
preview(vid); start(vid);
 
rgbImage = getdata(vid);
h=imshow(rgbImage,'Parent', handles.axes1);
% set(gcf, 'Position', get(0,'Screensize'));
fullImageFileName = fullfile(pwd, 'myfirstimage.jpg');
imwrite(rgbImage,fullImageFileName);
stop(vid); delete(vid);
[row col i]= size(rgbImage);
z=strcat(num2str(row),'x',num2str(col));
% h= imshow(img);
im = imagemodel(h);
q=  getImageType(im);

set(handles.edit2,'String','myfirstimage.jpg');
set(handles.edit3,'String',z);
set(handles.edit4,'String',q);
img=rgbImage;
I=img;
imshow(I,'Parent', handles.axes1);
handles.Pic = I;

% Choose default command line output for Example_GUI
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
       {'*.png;*.jpg;*.gif;*.tif', 'All Image Files (*.png, *.jpg, *.gif, *.tif)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Pick a file');        
            
display('loading....');display(filename);
I = fullfile(pathname, filename);
img = imread(I);
h=imshow(img,'Parent', handles.axes1);
[row col i]= size(img);
z=strcat(num2str(row),'x',num2str(col));
% h= imshow(img);
im = imagemodel(h);
q=  getImageType(im);

set(handles.edit2,'String',filename);
set(handles.edit3,'String',z);
set(handles.edit4,'String',q);
% img=im2double(img);

if (i==3) 
Img=rgb2gray(img);
imshow(I,'Parent', handles.axes1);
% handles.Pic = I;
handles.Pic = img;
end

% Choose default command line output for Example_GUI
        handles.output = hObject;

% Update handles structure
        guidata(hObject, handles);


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
% Choose default command line output for Example_GUI


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
% featmatch_ratio = str2num(handles.popupmenu4);
% handles.featmatch_ratio = featmatch_ratio;


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
