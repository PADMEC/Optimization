function [G,dG,vargout]=StructLimitState(x,varargin)

global Base_G Base_dG
global G_const linStress linStressComputed HYPLAS_flag

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

%input of the nonlin analysis
%vargs={area,props,fext,glb,modos,RBtype,POD_Z};
% ang  = props(1,:);
% comp = props(2,:);
% els  = props(3,:);

if nargout>1
    % Sensibility
    sens=1;
else
    sens=0;
    varargin{5}=abs(varargin{5});
end
if ~HYPLAS_flag
    %input of the nonlin analysis
    %varargin={area,props,fext,glb,modos,RBtype,POD_Z};
    if sens
        % -modos
        varargin{5}=-abs(varargin{5});
    else
        varargin{5}=abs(varargin{5});
    end
end

LinLoad = 1e-3;
[np,nv]=size(x);
%Load linear pre-computations
if np==0
    %Linear analysis
    x=eye(nv)*LinLoad;np=nv;
    %nb=length(varargin{1});
    %sig=zeros(nb,nv);
    precomputations=1;
else
    %sig=zeros(1,length(G_const));
    precomputations=0;
end
vargout={};

if (np>nv)
    x=x';
    [np,nv]=size(x);
end
    
linStressComputed=0;
if linStressComputed
    %% Linear Stress Base
    A=zeros(length(linStress{1}),nv);
    for iv=1:nv
        A(:,iv)=linStress{iv};
    end
    
    %dsl=zeros(length(G_const),nv);
    for ip=1:np
        
        % Linear Stress Comput
        LinSigs=A*x(ip,:)';
        test_lin(ip) = max(abs(LinSigs));
        
        % LimitStateFunc Base
        for ig=length(G_const)
            %% Linear Failure Function Comput
            B{ig}=zeros(length(Base_G{1,ig}),nv);
            for iv=1:nv
                B{ig}(:,iv)=Base_G{iv,ig};
            end
            C{ig}=B{ig}*x(ip,:)';
            
            %% Grad of Linear Failure Function Comput
            if var_type==1
                dC{ig}=Base_dG{1,ig}*x(ip,1);
                for iv=2:nv
                   dC{ig}=dC{ig} + Base_dG{iv,ig}*x(ip,iv);
                end
            else
                dC{ig}=B{ig};
            end

        end
        %% Failure Function Value
        [G,dG]=Lim_State_comput(C,dC,sens,precomputations,ip);
        
    end
else
    test_lin=1;
end

% global mat_curv
% test_nolin = max(test_lin)>mat_curv(1,2);
test_nolin = 1;
if test_nolin

    varargin0=varargin;
    for ip=1:np

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Run Analysis, G(U) computation
        if HYPLAS_flag
            [v_out]=Stoch_HyplasAnalysis(x(ip,:),varargin0{:});
        else
            [v_out]=Stoch_Analysis(x(ip,:),varargin0{:});
        end
        if sens
            nout=length(v_out)/2;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % TEST
%             ndv=length(x);
%             for i=1:ndv
%                 dx = zeros(1,ndv);
%                 dx(i)=x(i)/10000;
%                 [v_out1]=Stoch_Analysis(x(ip,:)+dx,varargin0{:});
%                 for iout=1:4
%                     Grads{iout}(i,:)=[(v_out1{iout}(1)-v_out{iout}(1))/dx(i)-v_out{iout+nout}(i,1),v_out{iout+nout}(i,1)];
%                 end
%             end
%             disp('StructLimitState')
%             Grads{:}
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            nout=length(v_out);
        end
        vargout=v_out(:);

        %G_const = 1 VOLUME, 2 DESLOCAMENTO, 3 TENSAO, 4 COMPILANCE,
        %           5 lambd
        %vol,u,vonmis,sn,lambd

        %Load linear pre-computations
        if precomputations
            linStress{ip}=Sigs(:,end)/LinLoad;
            linStressComputed=1;
        end
        ipg=0;
        for ig=G_const
            if ig>nout
                warning(['StructLimitState: G_const > num of output'...
                     'arguments of Stoch_Analysis'])
                 return
            end
            
            ipg=ipg+1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Failure Function type
            dC=[];
            C{ipg}=v_out{ig}; 
            if sens
                dC{ipg}=v_out{ig+nout}';
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if precomputations
                Base_G{ip,ipg} = C{ipg}/LinLoad;
                if sens
                    Base_dG{ip,ipg} = dC{ipg}/LinLoad;
                end
                G=0;dG=0;
            else
                [G,dG]=Lim_State_comput(C,dC,sens);
            end

        end

    end
end
    
function [G,dG,test_lin]=Lim_State_comput(C,dC,sens)

global G_const G_limits set_Failure_FORM ABS_Failure_func
ipg = 0;dG =0;
test_lin=0;
ABS_Failure_func = 0;%Regular formulation
%ABS_Failure_func = 1;%Paper error formulation case
for ig=G_const
    ipg=ipg+1;
    
    if isempty(set_Failure_FORM)
        if sign(G_limits)==-1
            [g(ipg),k_LS]=max((C{ipg}));
        else
            if ABS_Failure_func
                [g(ipg),k_LS]=max((C{ipg}));
            else
                [g(ipg),k_LS]=max(abs(C{ipg}));
            end
        end
    elseif set_Failure_FORM==0
        if (sign(G_limits)==-1)||(ABS_Failure_func)
            g(:,ipg)=C{ipg};
        else
            g(:,ipg)=abs(C{ipg});
        end
    else
        if (sign(G_limits)==-1)||(ABS_Failure_func)
            g(ipg)=C{ipg}(set_Failure_FORM);
        else
            g(ipg)=abs(C{ipg}(set_Failure_FORM));
        end
    end

    if sens
        if isempty(set_Failure_FORM)
            % Worst failure func
            if (sign(G_limits)==-1)||(ABS_Failure_func)
                dg(:,:,ipg) = dC{ipg}(k_LS,:);
            else
                sgn = sign(C{ipg}(k_LS));
                oneM=dC{ipg}(1,:)*0+1;
                dg(:,:,ipg) = (sgn*oneM).*dC{ipg}(k_LS,:);
            end
            
        elseif set_Failure_FORM==0
            % All failure funcs
            if (sign(G_limits)==-1)||(ABS_Failure_func)
                dg(:,:,ipg) = dC{ipg};
            else
                sgn = sign(C{ipg});
                oneM=dC{ipg}(1,:)*0+1;
                dg(:,:,ipg) = (sgn*oneM).*dC{ipg};
            end
        else
            % Specific failure func
            if (sign(G_limits)==-1)||(ABS_Failure_func)
                dg(:,:,ipg) = dC{ipg}(set_Failure_FORM,:);
            else
                sgn = sign(C{ipg}(set_Failure_FORM));
                sgn = sign(C{ipg}(set_Failure_FORM));
                oneM=dC{ipg}(1,:)*0+1;
                dg(:,:,ipg) = (sgn*oneM).*dC{ipg}(set_Failure_FORM,:);
            end
        end
            
    end
end

% g < G_limits (-)
% G = g-G_limits < (-)

G = g./abs(G_limits) - sign(G_limits);
%[G(ip),imax]=max(Gi);
if sens
    %dG = ds(imax,:)'/G_limits(imax);
    dG = dg'/abs(G_limits);
end