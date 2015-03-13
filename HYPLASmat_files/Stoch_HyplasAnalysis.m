function [v_out]=Stoch_HyplasAnalysis(x,varargin)
% [v_out]=Stoch_Analysis(x,varargin)
%
% Comput NL-FEM Analsysis outputs (and sensitives) for random variable x
% varargin_i=RandomStructDef(x,varargin{:});
% v_out = {vol,u,vonmis,sn,dvol,du,dsig,den};
% varargin={x, hyplasprojname,0,RandVarDist(2,:)};

global var_type sens

funame = 'Hyplas_analysis';

if isempty(sens)
    sens = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables Definitions (U)
varargin{4}=x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Analysis
[Fs,Us,Sigs,sn,vol]=feval(funame, varargin{:});
u = Us(:,end);
vonmis = Sigs(:,end);
v_out = {vol,u,vonmis,sn};

if sens
    if (var_type==1)
        % Design Variables
        id = 1;
    else
        % Random Variables
        id =4;
    end
    var0=varargin{id};
    [np,nv]=size(var0);
    nvout=length(v_out);
    for iv=1:nv
        var=var0;
        dv=var(iv)*1e-6;
        dv(dv==0)=1e-6;
        var(iv)=var(iv)+dv;
        varargin{id}=var;
        [Fs,Us,Sigs,sn,vol]=feval(funame, varargin{:});
%         u = Us(:,end);
%         vonmis = Sigs(:,end);
        v_out1 = {vol,Us(:,end),Sigs(:,end),sn};
        for ivout=1:nvout
            dv_out{ivout}(iv,:)=(v_out1{ivout}'-v_out{ivout}')/dv;
        end
    end
    v_out = {v_out{:},dv_out{:}};
end

    