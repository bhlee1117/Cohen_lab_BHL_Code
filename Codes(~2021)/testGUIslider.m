function varargout = testGUIslider(varargin)
% TESTGUISLIDER M-file for testGUIslider.fig
%      TESTGUISLIDER, by itself, creates a new TESTGUISLIDER or raises the existing
%      singleton*.
%
%      H = TESTGUISLIDER returns the handle to a new TESTGUISLIDER or the handle to
%      the existing singleton*.
%
%      TESTGUISLIDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTGUISLIDER.M with the given input arguments.
%
%      TESTGUISLIDER('Property','Value',...) creates a new TESTGUISLIDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testGUIslider_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testGUIslider_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testGUIslider

% Last Modified by GUIDE v2.5 21-Apr-2010 14:40:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testGUIslider_OpeningFcn, ...
                   'gui_OutputFcn',  @testGUIslider_OutputFcn, ...
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

% --- Executes just before testGUIslider is made visible.
function testGUIslider_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testGUIslider (see VARARGIN)

% Choose default command line output for testGUIslider
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes testGUIslider wait for user response (see UIRESUME)
% uiwait(handles.figure1);

t = 0:0.01:10*pi;
y = 5*sin(t);
plot(t, y)
grid on
axis([0 5*pi -5 5])



% --- Outputs from this function are returned to the command line.
function varargout = testGUIslider_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

temp1 = get(handles.slider1, 'Value');
set(handles.axes1, 'ylim', [-(5+temp1) 5+temp1])


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

temp2 = get(handles.slider2, 'Value');
set(handles.axes1, 'xlim', [0 (5+temp2)*pi])


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
