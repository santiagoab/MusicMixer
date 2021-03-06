function varargout = GUIv1(varargin)
% GUIV1 MATLAB code for GUIv1.fig
%      GUIV1, by itself, creates a new GUIV1 or raises the existing
%      singleton*.
%
%      H = GUIV1 returns the handle to a new GUIV1 or the handle to
%      the existing singleton*.
%
%      GUIV1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIV1.M with the given input arguments.
%
%      GUIV1('Property','Value',...) creates a new GUIV1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIv1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIv1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIv1

% Last Modified by GUIDE v2.5 07-Jun-2018 03:00:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIv1_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIv1_OutputFcn, ...
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


% --- Executes just before GUIv1 is made visible.
function GUIv1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIv1 (see VARARGIN)

% Choose default command line output for GUIv1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIv1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIv1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browseButton.
function browseButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[songs, folder] = getMusicFiles();
set(handles.dirText, 'String', folder);
assignin('base','songs',songs) %save var to base workspace
assignin('base','folder',folder)
features = doToSongs(folder,songs);
assignin('base','features',features)
sortedByEnergy = sortFeaturesEnergy(features);
assignin('base','sortedByEnergy',sortedByEnergy)
sortedByBPM = sortFeaturesBPM(features);
assignin('base','sortedByBPM',sortedByBPM)



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes on button press in Low.
function Low_Callback(hObject, eventdata, handles)
% hObject    handle to Low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Low


% --- Executes on button press in Medium.
function Medium_Callback(hObject, eventdata, handles)
% hObject    handle to Medium (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Medium


% --- Executes on button press in High.
function High_Callback(hObject, eventdata, handles)
% hObject    handle to High (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of High


% --- Executes on button press in Create.
function Create_Callback(hObject, eventdata, handles)
% hObject    handle to Create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
energyBtn = get(handles.buttongroup,'SelectedObject');
btnPressed = get(energyBtn,'String');
assignin('base','btnPressed',btnPressed)
sortedByEnergy = evalin('base', 'sortedByEnergy');
playlistSongs = createPlaylist(sortedByEnergy,btnPressed);
assignin('base','playlistSongs',playlistSongs)
