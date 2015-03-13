function [Mu,Cov]=check_DesRand(Mu,Cov,varargin)
% Design and Random Variable Check!

global RandVarProp rv_data_comp link

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

[~,nv]=size(Mu);
chang_count=0;
cof = Cov./(Mu'*Mu);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Variables Definitions (U)
    for iv=1:nv
        RVkind = RandVarProp(1,iv);
        RVwich = RandVarProp(2,iv);
        if (RVkind == 1)&&(RVwich == 1)
            % Design and Random Variable!
            chang_count=1;
            components=rv_data_comp{iv};
            ic=find(link==components(1),1);
            Mu(iv) = varargin{RVkind}(ic);
        end
    end
    if chang_count
        Cov = cof.*((Mu'*Mu));
    end
end