function [f,df]=constfun(x,funame,~,inv,~,~,l,~,varargin)

global GradOutput_on
if GradOutput_on
    [f,df]=feval(funame,x,varargin{:});
    M_l = ones(length(x),1)*l(inv);
    df=sum( (df(:,inv))./M_l ,2);
else
    f=feval(funame,x,varargin{:});
end
l(l==0)=1;l=abs(l);
f=sum(f(inv)./l(inv));
