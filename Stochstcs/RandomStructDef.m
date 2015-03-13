function varargin=RandomStructDef(x,varargin)

global RandVarProp rv_data_comp 

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

[ip,nv]=size(x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Variables Definitions (U)
    for iv=1:nv
        components=rv_data_comp{iv};
        RVkind = RandVarProp(1,iv);
        RVwich = RandVarProp(2,iv);
        updated_varg = RV_update(x(ip,iv),RVkind,RVwich,...
            components,varargin{RVkind});
        varargin{RVkind} = updated_varg;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %varagout=varargin;
end