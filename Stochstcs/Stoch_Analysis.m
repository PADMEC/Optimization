function [v_out]=Stoch_Analysis(x,varargin)
% [v_out]=Stoch_Analysis(x,varargin)
%
% Comput NL-FEM Analsysis outputs (and sensitives) for random variable x
% varargin_i=RandomStructDef(x,varargin{:});
% v_out=NLfesol(varargin_i{:});
% v_out = {vol,u,vonmis,sn,dvol,du,dsig,den};


global linStress linStressComputed

sens = 1;
LinLoad = 1e-3;
[np,nv]=size(x);
%Load linear pre-computations
if np==0
    %Linear analysis
    x=eye(nv)*LinLoad;np=nv;
    precomputations=1;
else
    precomputations=0;
end

%% Linear Stress Base
linStressComputed=0;
if linStressComputed
    A=zeros(length(linStress{1}),nv);
    for iv=1:nv
        A(:,iv)=linStress{iv};
    end
    
    %dsl=zeros(length(G_const),nv);
    % Linear Stress Comput
    LinSigs=A*x(ip,:)';
    test_lin(ip) = max(abs(LinSigs));

    % LimitStateFunc Base
    %% Linear Failure Function Comput
    B{ig}=zeros(length(Base_G{1}),nv);
    for iv=1:nv
        B{ig}(:,iv)=Base_G{iv};
    end
    C{ig}=B{ig}*x(ip,:)';

    %% Grad of Linear Failure Function Comput
    if var_type==1
        dC=Base_dG{1}*x(ip,1);
        for iv=2:nv
           dC=dC + Base_dG{iv}*x(ip,iv);
        end
    else
        dC=B;
    end
else
    test_lin=inf;
end
% global mat_curv 
% testlin=max(test_lin)>mat_curv(1,2);
    testlin=1;
if testlin

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Variables Definitions (U)
    varargin_i=RandomStructDef(x,varargin{:});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run Analysis
    if sens
        varargin_i{5}=-abs(varargin{5});
        [Fs,Us,Sigs,sn,vol,T,floc,u_flamb,Lambd,D,freq,dF,du,dsig,den,dvol,dLamb]=NLfesol(varargin_i{:});
        %vargout={Fs,Us,Sigs,en,vol,T,floc,u_flamb,Lambd,D,freq,dF,du,dsig,den,dvol};
    else
        varargin_i{5}=abs(varargin{5});
        [Fs,Us,Sigs,sn,vol,T,floc,Lambd]=NLfesol(varargin_i{:});
        %vargout={Fs,Us,Sigs,en,vol,T,floc};
    end

    u = Us(:,end);
    vonmis = Sigs(:,end);
    v_out = {vol,u,vonmis,sn,Lambd};
    if sens, v_out = {v_out{:},dvol,du',dsig',den',dLamb};end
    

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % TEST
%     global var_type link
%     if var_type==2 % X = Rand var
%     else% X = Design var
%         varargin=varargin_i;
%         area = varargin{1};
%         ndv=max(link(:,1));
%         x = zeros(1,ndv);
%         for ivv=1:ndv
%             kv=find(link(:,1)==ivv,1);
%             x(ivv)= area(kv);
%         end
%     end
%     ndv=length(x);
%     for i=1:ndv
%         dx = x*0;
%         dx(i)=x(i)/10000;
%         dx(x==0)=1e-6;
%         if var_type==2 % X = Rand var
%         varargin_i=RandomStructDef(x+dx,varargin{:});varargin_i{5}=-abs(varargin{5});
%         else % X = Design var
%             varargin_i=varargin;
%             varargin_i{5}=-abs(varargin{5});
%             varargin_i{1}(link(:,1)==i)=x(i)+dx(i);
%         end
%         [Fs,Us,Sigs,sn,vol,T,floc,Lambd]=NLfesol(varargin_i{:});
%         v_out1 = {vol,Us(:,end),Sigs(:,end),sn,Lambd};
%         for iout=1:4
%             Grads{iout}(i,:)=[(v_out1{iout}(1)-v_out{iout}(1))/dx(i),v_out{iout+5}(i,1)];
%         end
%     end
%     disp('Stoch_Analysis')
%     Grads{:}
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Load linear pre-computations
    if precomputations
        linStress{ip}=Sigs(:,end)/LinLoad;
        linStressComputed=1;
    end
end