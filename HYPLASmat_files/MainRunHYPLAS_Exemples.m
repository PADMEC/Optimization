
clear variables global


%path('bin',path)
path('HYPLASgui_matlab',path)

callHYPLAS = 'bin\hyplas90a.exe';


%x = [2, 4.5];

% hObject    handle to pushbutton_visualise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%
% Variable shared with HYPLASguiPRE and HYPLASguiPOST
global hyplasprojname
% Variables shared with PLOT RESULTS pushbutton funtion of HYPLASguiPOST
global lnods coord incrtable availflag loadfac filenameres...
       displ reac  nodepresc varnames  varvalues ngrup

% Retrieve HYPLAS results file name
%hyplasprojname='13_6_2';
%x = [4,4]
hyplasprojname='KubrickSqrtHole';
x = [250,250]


%%% Exemples Dat
exemplesdir='book_examples/data_files';
cd(exemplesdir)
datfiles=dir('*.dat');
cd('..\..')
ndt = length(datfiles);

for i=1:ndt
    hyplasprojname = datfiles(i).name(1:end-4);

    name_res=[ hyplasprojname '.res']; 
    filenameres=[exemplesdir '/' name_res];
    name_dat=[ hyplasprojname '.dat']; 
    filenamedat=[exemplesdir '/' name_dat]

    str_command = sprintf('%s %s 0',callHYPLAS,filenamedat);

    system(str_command,'-echo')

    exres=exist(filenameres,'file');
    try
    if exres
        [lnods,    coord,     incrtable, availflag,...
         loadfac,  displ,     reac,      nodepresc,...
         varnames, varvalues, ngrup                  ] = hyplaspost(filenameres);
    else
        filenamedat=[hyplasprojname '.dat'];exdat=exist(filenamedat,'file');
        [lnods,    coord,     incrtable, availflag,...
         loadfac,  displ,     reac,      nodepresc,...
         varnames, varvalues, ngrup                  ] = hyplaspost(filenameres);
    end

    %
    %1 - reactions, 2 - contours, 3 - Deformed mesh, 4 - Initial mesh, 5 - displacements,
    % 6 - incremental contours, 7 - incremental displacements, 8 - incremental reactions, 9,
    %10, 11 - Node numbers
    checkbox=zeros(1,11);
    checkbox([1,2,3])=1;
    PlotResults(checkbox);
    catch
        disp('errorrrr!!!!') 
        disp(filenameres)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%