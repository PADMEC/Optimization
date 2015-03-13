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
function [lnods,    coord,     incrtable, availflag,...
          loadfac,  displ,     reac,      nodepresc,...
          varnames, varvalues, ngrup ,totarea] = hyplaspost(filenameres)
%
%
%
% Open HYPLAS results output file
fid = fopen(filenameres, 'r');
%
% Check that the file has only one group of elements
for iloop=1:1000000
    tline = fgetl(fid);
    if length(tline)>22 && strcmp(tline(2:15),'Element Groups')
        break
    end
end
% Total number of element groups
ngrup = str2double(tline(52:length(tline)));
%
% Find and read connectivity table
% --------------------------------
for iloop=1:1000000
    tline = fgetl(fid);
    if length(tline)>22 && strcmp(tline(2:23),'Element connectivities')
        break
    end
end
% read total nunber of elements
nelem = str2num(tline(52:length(tline)));
% read connectivity table
for iloop=1:4
    tline = fgetl(fid);
end
%
tline = fgetl(fid);
lnods(:,1) = str2num(tline);
nnode = length(lnods)-2;
%
if (nelem > 1)
    if nnode == 3
        lnods(:,2:nelem) = fscanf(fid, '\n%d %d %d %d %d', [nnode+2,nelem]);
    elseif nnode == 4
        lnods(:,2:nelem) = fscanf(fid, '\n%d %d %d %d %d %d', [nnode+2,nelem]);
    elseif nnode == 8
        lnods(:,2:nelem) = fscanf(fid, '\n%d %d %d %d %d %d %d %d %d %d', [nnode+2,nelem]);
    else
        error('error in hyplaspost - element type not implemented')
    end
else
    % Do a dummy scan
    fscanf(fid, '\n%d', [1,1]);
end

lnods = lnods'; elorder = lnods(:,1);
lnods = lnods(:,3:nnode+2); lnods(elorder,:) = lnods;
%
% Find and read initial nodal coordinates
% ---------------------------------------
for iloop=1:100000
    tline = fgetl(fid);
    if length(tline)>23 && strcmp(tline(1:24),'Nodal point co-ordinates')
        break
    end
end
% read total nunber of nodes
npoin = str2num(tline(52:length(tline)));
% read coordinates
for iloop=1:4
    tline = fgetl(fid);
end
coord.init = fscanf(fid, '\n%f %f %f', [3, npoin]);
coord.init = coord.init'; nodorder = coord.init(:,1);
coord.init = coord.init(:,2:3); coord.init(nodorder,:) = coord.init;
%
% number of nodes with prescribed displacements
for iloop=1:100000
    tline = fgetl(fid);
    if length(tline)>25 && strcmp(tline(1:25),' Prescribed displacements')
        break
    end
end
npresc = str2num(tline(77:length(tline)));

%
% initialize reactions arrays
reac.x = zeros(npresc,1); reac.y = zeros(npresc,1); nodepresc = zeros(npresc,1);
reac.all.x = zeros(npoin,1);reac.all.y=zeros(npoin,1);
% and tables of increments available for plotting, corresponding load
% factors and variable availability flags
incrtable = []; loadfac = []; availflag = []; varvalues = [];
%
% Read most of plotting data
%
firsttimeflag=0;varnames={};
while ~(isnumeric(tline) && tline==-1) % loops until it hits end of file
    tline = fgetl(fid);
    tlleng=length(tline);
    % read increment information
    if tlleng>26 && strcmp(tline(1:22),' Results for increment');
        iincr = str2num(tline(24:27)); lf = str2num(tline(63:tlleng));
        loadfac = [loadfac lf]; incrtable = [ incrtable iincr ];
    % read displacements
    elseif tlleng>25 && strcmp(tline(1:26),' Displacement of structure')
        availflag(1,iincr) = 1; % set displacements availability flag
        for ii=1:3
            tline = fgetl(fid);
        end
        displaux = fscanf(fid, '\n%f %f %f', [3, npoin]);
        displaux = displaux'; nodorder = displaux(:,1);
        displaux = displaux(:,2:3); displaux(nodorder,:) = displaux;
        displ.x(:,iincr) = displaux(:,1); displ.y(:,iincr) = displaux(:,2);
    
    elseif length(tline)>14 && strcmp(tline(1:14),'  Total Area =')
        totarea=sscanf(tline,'  Total Area = %g');

    % read reactions
    elseif tlleng>9 && strcmp(tline(1:10),' Reactions') && ~( tlleng>11 && strcmp(tline(1:12),' Reactions..') );
        availflag(2,iincr) = 1; % set reactions availability flag
        for iloop=1:3
            tline = fgetl(fid);
        end
        [~,nstring]=sscanf(tline,'%s ');
        displaux = fscanf(fid, '\n%d %f %f', [nstring,npresc])';
        inode = displaux(:,1); nodepresc(:,iincr) = inode;
        reac.x(:,iincr)=displaux(:,2); reac.y(:,iincr)=displaux(:,3);
        reac.all.x(inode,iincr) = displaux(:,2);
        reac.all.y(inode,iincr) = displaux(:,3);
%         for ipresc = 1:npresc
%             tline = fgetl(fid);
%             inode = str2num(tline(1:5)); nodepresc(ipresc,iincr) = inode;
%             xforce = str2num(tline(6:22)); yforce = str2num(tline(23:39));
%             reac1.x(ipresc,iincr) = xforce; reac1.y(ipresc,iincr) = yforce;
%             reac1.all.x(inode,iincr) = xforce; reac1.all.y(inode,iincr) = yforce;
%         end
        
    % Read averaged nodal variables for countour plot.
    % Store variable names in varnames and all variables values for all
    % nodes and all increments in varvalues
    elseif tlleng>20 && strcmp(tline(2:21),'Gauss point stresses') && ngrup==1
        get_gauss
    elseif tlleng>14 && strcmp(tline(1:15),' Averaged nodal') && ngrup==1
        get_averages
    end
end

status = fclose(fid);