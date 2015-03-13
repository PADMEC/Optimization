clear variables global
% Modified by Renato S. Motta (03/2014)
global bindir HYPLASexe
global NO_IG
NO_IG=1;
%path('bin',path)
path('HYPLASgui_matlab',path)
path('Mfiles-NonL-POD',path)
bindir = 'bin';
HYPLASexe = 'hyplas90a.exe';



% hObject    handle to pushbutton_visualise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%
% Variable shared with HYPLASguiPRE and HYPLASguiPOST
global hyplasprojname
% Variables shared with PLOT RESULTS pushbutton funtion of HYPLASguiPOST
global lnods coord incrtable availflag loadfac ...
       displ reac  nodepresc varnames  varvalues ngrup

% Retrieve HYPLAS results file name
hyplasprojname='13_6_2';

x = [4, 4];
[Fs,Us,Sigs,ien,vol] = Hyplas_analysis( x, hyplasprojname);

[Fs2,Us2,Sigs2,ien2,vol2,...
    exres,lnods,    coord,     incrtable, availflag,...
     loadfac,  displ,     reac,      nodepresc,...
     varnames, varvalues, ngrup     ] = Hyplas_analysis( x, hyplasprojname);
 
%
%1 - reactions, 2 - contours, 3 - Deformed mesh, 4 - Initial mesh, 5 - displacements,
% 6 - incremental contours, 7 - incremental displacements, 8 - incremental reactions, 9,
%10, 11 - Node numbers
checkbox=zeros(1,11);
checkbox([2,4])=1;
PlotResults(checkbox);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%