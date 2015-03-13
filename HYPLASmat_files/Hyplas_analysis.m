function [varargout]=Hyplas_analysis(x, hyplasprojname1,RBM,randvar)
%Fs,Us,Sigs,arq.ien,arq.ivol,T
%RBM==-3 for FEM using command input
%RBM== 0 for FEM analysis
%RBM > 0 for RBM analysis (1 - POD, > - classic RBM)
%RBM < 0 for analysis (FEM or RBM) with command string given (hyplasprojname1)

global bindir HYPLASexe   Stochast_on
% hyplasprojname
% command_string

global gambi
if gambi==1
    randvar(3)=x(3);
    x(3)=[];
end
if nargin>2
    % RBM parameter data given..
    if RBM>0, str_num_args = sprintf(' %.12g ',[-length(x) x]);
    elseif (RBM==-3)||(RBM==0)
        % FEM using command input
        str_num_args = sprintf(' %.12g ',[length(x) x]);
    else
        % -1 or -2
        str_num_args = sprintf(' %.12g ',[-length(x) x]);
    end
    
    if Stochast_on
        str_num_args = [str_num_args,sprintf(' %.12g ',[length(randvar) randvar])];
    end
    
    if RBM<0
        % Using hyplasprojname1 as command_string input data
        command_string=hyplasprojname1;
        if isempty(command_string)
            disp('Error.. "command_string" not defined!')
            pause()
            return
        end
        str_command = [command_string ,' ',str_num_args];
        system(str_command,'-echo');
        filenameres = 'HyplasResults.bin';
        [Fs,Us,Sigs,energy,vol] = GetHYPLASres(filenameres);
        varargout = {Fs,Us,Sigs,energy(end),vol};
        return
    end
else
    str_num_args = sprintf(' %.12g ',[length(x) x]);
    if Stochast_on
        str_num_args = [str_num_args,sprintf(' %.12g ',[length(randvar) randvar])];
    end
end

hyplasprojname = hyplasprojname1;
callHYPLAS = [bindir '/' HYPLASexe];


%     bindir1=strrep(bindir,'C:','/cygdrive/c');
%     bindir1=strrep(bindir1,'\','/');

name_res=[ hyplasprojname '.res']; filenameres=[ name_res];
name_dat=[ hyplasprojname '.dat']; filenamedat=[ name_dat];

if nargout<6
    %% Post-Proc for optimization compatibility
    command_string = sprintf('%s %s -resume_off -echo_off -debug_off ',callHYPLAS,filenamedat);
    if nargout==1
        % Asking just for the command_string
        varargout = {command_string};
        return
    end
    str_command = [command_string ,' ',str_num_args];
    system(str_command,'-echo');
    
    filenameres = 'HyplasResults.bin';
    [Fs,Us,Sigs,energy,vol] = GetHYPLASres(filenameres);
    varargout = {Fs,Us,Sigs,energy(end),vol};
    return
    
else
    %% Post-Proc for plot compatibility
    str_command = sprintf('%s %s -resume_on -echo_on -debug_off %s',callHYPLAS,filenamedat,str_num_args);
    system(str_command,'-echo');

    exres=exist(filenameres,'file');
    if exres
        [lnods,    coord,     incrtable, availflag,...
         loadfac,  displ,     reac,      nodepresc,...
         varnames, varvalues, ngrup , totarea ] = hyplaspost(filenameres);
    else
        exdat=exist(filenamedat,'file');
        if exdat
            [lnods,    coord,     incrtable, availflag,...
             loadfac,  displ,     reac,      nodepresc,...
             varnames, varvalues, ngrup , totarea ] = hyplaspost(filenameres);
        else
            disp(['Files not found!' filenameres])
        end
    end
    varargout = {exres, lnods,    coord,     incrtable, availflag,...
       loadfac,  displ,     reac,      nodepresc,...
     varnames, varvalues, ngrup};

% %% Post-Proc to optimization compatibility
%  incrid = incrtable(end);
%  %Fs,Us,Sigs,arq.ien,arq.ivol,T
%  nn  = length(displ.x);
%  Fs = zeros(2*nn,incrid);Us = Fs; 
%  %Sigs = zeros(nn,incrid);
%  Fs(1:2:end,:) = reac.all.x;
%  Fs(2:2:end,:) = reac.all.y;
%  Us(1:2:end,:) = displ.x;
%  Us(2:2:end,:) = displ.y;
%  Sigs = varvalues(:,8,:); % S-eff
%  
%  % External Energy X = displ.x(nodepresc(:,1),2)'*reac.x(:,2)
%  % External Energy Y = displ.y(nodepresc(:,1),2)'*reac.y(:,2)
%  energy = Fs(:,incrid)'*Us(:,incrid);
%  vol = totarea;
%  varargout = {Fs,Us,Sigs,energy,vol};
end