%**********************************************************************
%
% Function HYPLASpost: A HYPLAS post-processor for HYPLASgui
% 
%
% DESCRIPTION
%
% This function reads a HYPLAS results file (output of HYPLAS v2.0/v2.0.1)
% and organises the relevant data in arrays, etc, suitable for plotting by
% the graphical interface program HYPLASgui.
%
%
% NOTE
%
% This version only reads files with one single group of elements.
%
%
%
% HISTORY
% EA de Souza Neto, April 2011:   Initial coding
%
%***********************************************************************
%
function [F,dis,Stress,energ,vol] = GetHYPLASres(filenameres)
%
%
%
%Check!
exres=exist(filenameres,'file');
if ~exist(filenameres,'file')
    disp(['File ' filenameres ' dos not exist (cant find it in GetHYPLASres)'])
    return
end
    
% Open HYPLAS results output file
fid = fopen(filenameres, 'r');
% Get Number of increments (Arc_length = 0)
nincrm=sscanf(fgetl(fid),'%d ');
% Get Volume
vol=[];
while isempty(vol)
    vol=sscanf(fgetl(fid),'%f ');
end

if nincrm>0
    F=zeros(1,nincrm);iter=zeros(1,nincrm);
end
i=0;
tline=fgetl(fid);
while tline~=-1
    i=i+1;
    if ((i==1)&&(nincrm>0))
        initi=1;
    else
        initi=0;
    end
    %Iteration, time, Fincr, fcount to converge
	jumpline=1;
    while jumpline
        a=sscanf(tline,'%d %f %f %d ');
        if isempty(a)
            tline=fgetl(fid);
        else
            jumpline=0;
        end
    end
    iter(i)=a(1);
    F(i)=a(3);
    
    % Displacements
    ndis=sscanf(fgetl(fid),'%d ');
    if initi, dis=zeros(ndis,nincrm);end
    dis(1:ndis,i)=sscanf(fgetl(fid),'%f ');
    
    %Effective Stress - Gauss Point 
    sizes=sscanf(fgetl(fid),'%d '); % (n pont_gauss, n_elem) 
    ngp=sizes(1);
    if initi
        StressGP=zeros(sizes(2),ngp,nincrm);
        Stress=zeros(sizes(2)*ngp,nincrm);
    end
    for igp=1:ngp
        StressGP(:,igp,i)=sscanf(fgetl(fid),'%f ');
        Stress(igp:ngp:end,i)=StressGP(:,igp,i);
    end
    
    %Effective Stress - Nodal Average
    npoints=sscanf(fgetl(fid),'%d ');
    if initi, StressAv=zeros(npoints,nincrm);end
    StressAv(1:npoints,i)=sscanf(fgetl(fid),'%f ');
    
    if initi, energ=zeros(1,nincrm);end
    energ(i)=sscanf(fgetl(fid),'%f ');
    
    %Try new iteration
    tline=fgetl(fid);
end
fclose(fid);