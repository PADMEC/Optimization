function varargout = metodomo(varargin)
% METODOMO M-file for metodomo.fig
%      METODOMO, by itself, creates a new METODOMO or raises the existing
%      singleton*.
%
%      H = METODOMO returns the handle to a new METODOMO or the handle to
%      the existing singleton*.
%
%      METODOMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in METODOMO.M with the given input arguments.
%
%      METODOMO('Property','Value',...) creates a new METODOMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before metodomo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to metodomo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help metodomo

% Last Modified by GUIDE v2.5 22-May-2006 19:51:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @metodomo_OpeningFcn, ...
                   'gui_OutputFcn',  @metodomo_OutputFcn, ...
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


% --- Executes just before metodomo is made visible.
function metodomo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to metodomo (see VARARGIN)

% Choose default command line output for metodomo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes metodomo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = metodomo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
uiwait
h=guidata(gcf);
varargout{1} = h.fs;
            close;

% --- Executes on button press in sm.
function sm_Callback(hObject, eventdata, handles)
% hObject    handle to sm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
            handles.fs=1;
            guidata(gcf,handles);
            uiresume;

% --- Executes on button press in mm.
function mm_Callback(hObject, eventdata, handles)
% hObject    handle to mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 'Min-Max',
            handles.fs=2;
            guidata(gcf,handles);
            uiresume;

% --- Executes on button press in nnc.
function nnc_Callback(hObject, eventdata, handles)
% hObject    handle to nnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 'NNC',
            handles.fs=4;
            guidata(gcf,handles);
            uiresume;

% --- Executes on button press in nbi.
function nbi_Callback(hObject, eventdata, handles)
% hObject    handle to nbi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 'NBI',
            handles.fs=3;
            guidata(gcf,handles);
            uiresume;

% --- Executes on button press in fatias.
function fatias_Callback(hObject, eventdata, handles)
% hObject    handle to fatias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 'Fatias',
            handles.fs=5;
            guidata(gcf,handles);
            uiresume;


