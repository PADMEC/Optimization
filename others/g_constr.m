
function C = g_constr(U,varargin)

    %10/x+10/y-a=0; 10/x=a-10/y, x=1/(a/10-1/y)
    %y=[0.5,1.5], x=[1/(5-)]
%     y=6.5:.1:29;x=1./(1/5-1./y);
%     hold on,plot(x,y), plot(Mu(1),Mu(2),'+k')

% Input random variable data
%     Mu=[13,15];
%     Su=[2.,3.];
%     C=[1, -.6;-.6, 1];
%     Cu = C.*(Su'*Su);
%     ft=[1,1];
%     %DinamcVar = 30000;Method='MC';
%     DinamcVar = Mu;Method='FORM';
%     G_fun = 'g_constr';
%     
    C = 10./U(:,1)*100+10./U(:,2) - 2;

end
