function [gV,dgV]=PMA_reliab_obj(V,Mu,Su,Jzu,ft,nrv,parmeters,G_fun,varargin)

% global betatarg

% V - standad space (RV)
% U - original space (RV)

    U = transNf(V,Jzu,ft,nrv,parmeters);
    [gV,dgU] = feval(G_fun,U(:)',varargin{:});
    Jzv = inv_normeq(Mu,Su,U,ft);
    dgV = Jzu*Jzv*dgU;
     gV =  -gV;
    dgV = -dgV;
%     
%     gV =  gV*1e5;
%     dgV= dgV*1e5;
end