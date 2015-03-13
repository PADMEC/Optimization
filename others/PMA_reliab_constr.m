function [C,d,dC,dd]=PMA_reliab_constr(V,Mu,Su,Jzu,ft,nrv,parmeters,G_fun,varargin)

% d=[];dd=[];
C=[];dC=[];

 global betatarg

% V - standad space (RV)
% U - original space (RV)
beta2=V'*V;
% C = [betatarg - beta; beta - betatarg];
d = [betatarg^2 - beta2];
% if beta==0;beta=1;end
% dC = [-2*V, 2*V];
dd = -2*V;

end