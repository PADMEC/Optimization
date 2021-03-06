function plot_HYPLAS(funame,varargin0)
%path('bin',path)
% path('HYPLASgui_matlab',path)

% callHYPLAS = 'bin\hyplas90a.exe';

%x = [2, 4.5];
%x = [3,3]

% hObject    handle to pushbutton_visualise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%
% Variable shared with HYPLASguiPRE and HYPLASguiPOST
% global hyplasprojname
% Variables shared with PLOT RESULTS pushbutton funtion of HYPLASguiPOST
global lnods coord incrtable availflag loadfac Sigs ...
       displ reac  nodepresc varnames  varvalues ngrup

% Retrieve HYPLAS results file name
%hyplasprojname='13_6_2';

% hyplasprojname='KubrickSqrtHole';
%x = [477,477]
% x = [200,200]

% tic
 [Fs,Us,Sigs,energy,vol]=feval(funame, varargin0{:});
% 
% toc

tic
[exres,lnods,    coord,     incrtable, availflag,...
     loadfac,  displ,     reac,      nodepresc,...
     varnames, varvalues, ngrup   ] = feval(funame, varargin0{:});
toc
if exres<0
    print(['Error... Output file not found: ' varargin0{2}])
    return
end
%
%1 - reactions, 2 - contours, 3 - Deformed mesh, 4 - Initial mesh, 5 - displacements,
% 6 - incremental contours, 7 - incremental displacements, 8 - incremental reactions, 9,
%10, 11 - Node numbers
checkbox=zeros(1,11);
% checkbox(4)=1;
checkbox([2,3])=1;
loadstep=0;
varchoice=0;
histplot=0;
if ~histplot
    PlotResults(checkbox,loadstep,varchoice);
else
    for i=incrtable
        PlotResults(checkbox,i,varchoice);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% name_res=[ hyplasprojname '.res']; filenameres=['bin/' name_res];
% name_dat=[ hyplasprojname '.dat']; filenamedat=['bin/' name_dat];
% 
% str_num_args = sprintf('%g ',[length(x) x]);
% str_command = sprintf('%s %s -resume_off -echo_off -debug_off %s',callHYPLAS,filenamedat,str_num_args)
% tic
% system(str_command,'-echo')
% toc
% 
% 
% filerests = 'HyplasResults.bin';
% [Fs,Us,Sigs,energy,vol] = GetHYPLASres(filerests);
%     
%     
% str_command = sprintf('%s %s -resume_on -echo_off -debug_on %s',callHYPLAS,filenamedat,str_num_args)
% tic
% system(str_command,'-echo')
% toc
% exres=exist(filenameres,'file');
% if exres
%     [lnods,    coord,     incrtable, availflag,...
%      loadfac,  displ,     reac,      nodepresc,...
%      varnames, varvalues, ngrup                  ] = hyplaspost(filenameres);
% else
%     filenamedat=[hyplasprojname '.dat'];exdat=exist(filenamedat,'file');
%     [lnods,    coord,     incrtable, availflag,...
%      loadfac,  displ,     reac,      nodepresc,...
%      varnames, varvalues, ngrup                  ] = hyplaspost(filenameres);
% end

