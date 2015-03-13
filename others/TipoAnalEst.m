function varargout = TipoAnalEst(varargin)
% TIPOANALEST M-file for TipoAnalEst.fig
%      TIPOANALEST, by itself, creates a new TIPOANALEST or raises the existing
%      singleton*.
%
%      H = TIPOANALEST returns the handle to a new TIPOANALEST or the handle to
%      the existing singleton*.
%
%      TIPOANALEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIPOANALEST.M with the given input arguments.
%
%      TIPOANALEST('Property','Value',...) creates a new TIPOANALEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TipoAnalEst_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TipoAnalEst_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help TipoAnalEst

% Last Modified by GUIDE v2.5 22-Sep-2006 17:30:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TipoAnalEst_OpeningFcn, ...
                   'gui_OutputFcn',  @TipoAnalEst_OutputFcn, ...
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


% --- Executes just before TipoAnalEst is made visible.
function TipoAnalEst_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TipoAnalEst (see VARARGIN)

% Choose default command line output for TipoAnalEst
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TipoAnalEst wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TipoAnalEst_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

uiwait
h=guidata(gcf);

h.modos=[1 h.modos].*h.tipo;

varargout{1} = h.modos;
            close;


% --- Executes on button press in flamb.
function flamb_Callback(hObject, eventdata, handles)
% hObject    handle to flamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flamb

v=get(hObject,'Value');
objd=findobj(gcf,'Tag','modflam');
if v==1
    set(objd,'enable','on')
else
    set(objd,'enable','off')
end

% --- Executes on button press in vib.
function vib_Callback(hObject, eventdata, handles)
% hObject    handle to vib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vib

v=get(hObject,'Value');
objd=findobj(gcf,'Tag','modvib');
if v==1
    set(objd,'enable','on')
else
    set(objd,'enable','off')
end


% --- Executes on button press in statc.
function statc_Callback(hObject, eventdata, handles)
% hObject    handle to statc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of statc




% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


htpa=findobj('tag','tipoa');
htps=get(htpa,'children');
nt=length(htps);
val=zeros(1,3);
modos=[0 0];
for i=1:nt
    name=get(htps(i),'tag');
    itp=nt;
    switch name
        case 'statc'
            itp=1;
        case 'flamb'
            itp=2;
            modf=ptn('modflam');
            modos(1)=modf;
        case 'vib'
            itp=3;
            modv=ptn('modvib');
            modos(2)=modv;
    end
    if get(htps(i),'value')==1
        val(itp)=get(htps(i),'value');
    end
end
%get(htpa)
            handles.tipo=val;
            handles.modos=modos;
            guidata(gcf,handles);
            uiresume;




function modflam_Callback(hObject, eventdata, handles)
% hObject    handle to modflam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modflam as text
%        str2double(get(hObject,'String')) returns contents of modflam as a double


% --- Executes during object creation, after setting all properties.
function modflam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modflam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function modvib_Callback(hObject, eventdata, handles)
% hObject    handle to modvib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modvib as text
%        str2double(get(hObject,'String')) returns contents of modvib as a double


% --- Executes during object creation, after setting all properties.
function modvib_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modvib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


