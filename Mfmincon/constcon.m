function [R,C,dR,D]=constcon(x,funame,coname,inv,fx,fm,l,par,varargin)

C=[];
D=[];

global GradOutput_on
if GradOutput_on
    [R,~,dR]=feval(coname,x,varargin{:});
    [f,df]=feval(funame,x,varargin{:});

    M_l = ones(length(x),1)*l(par);
    dres=df(:,par)./M_l;
    dR=[dR dres];
else
    [R]=feval(coname,x,varargin{:});
    [f]=feval(funame,x,varargin{:});
end
l(l==0)=1;l=abs(l);
res=(f(par)'-fx(:))./l(par)';
R=[R(:); res];