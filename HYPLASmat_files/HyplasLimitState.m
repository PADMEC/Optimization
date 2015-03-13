function [G,dG,vargout]=HyplasLimitState(x,varargin)

% varargin={x, hyplasprojname,0,RandVarDist(2,:)};
global G_const G_limits set_Failure_FORM

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

%input of the nonlin analysis
%vargs={area,props,fext,glb,modos,RBtype,POD_Z};
% ang  = props(1,:);
% comp = props(2,:);
% els  = props(3,:);

%input of the nonlin analysis
%varargin={area,props,fext,glb,modos,RBtype,POD_Z};
if nargout>1
    % Sensibility
    sens=1;
else
    sens=0;
end

[np,nv]=size(x);

vargout={};

varargin0=varargin;
for ip=1:np

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run Analysis, G(U) computation
    [vargout]=Stoch_HyplasAnalysis(x(ip,:),varargin{:});
    % vargout = {vol,u,vonmis,sn,dvol,du,dsig,den};

    %% Multiple output computation
    if iscell(vargout)
        n_out=length(vargout);
        n_grad = n_out/2;
        n_out = n_grad;
        v_out={vargout{1:n_out}};
        dv_out={vargout{n_out+(1:n_grad)}};
        %nvar=size(dv_out0{1},2);
    else
        disp('Output from Stoch_HyplasAnalysis not recognized')
    end
    %G_const = 1 VOLUME, 2 DESLOCAMENTO, 3 TENSAO, 4 COMPILANCE

    ipg=0;
    %% Compute Limite State
    for ig=G_const
        ipg=ipg+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Failure Function type
        if set_Failure_FORM>0
            kmax=set_Failure_FORM;
            C{ipg}=abs(v_out{ig}(kmax));
        else
            [C{ipg},kmax]=max(abs(v_out{ig}));
        end
        dC{ipg} = dv_out{ig}(:,kmax)*sign(v_out{ig}(kmax));
        G(ipg) = C{ipg}./G_limits - 1;
        %[G(ip),imax]=max(Gi);
        if sens
            %dG = ds(imax,:)'/G_limits(imax);
            dG(:,ipg) = dC{ipg}/G_limits;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

end
end