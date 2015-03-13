function [f,df]=wsfun(x,funame,coname,t,fm,l,a,n,par,varargin)
global GradOutput_on

if GradOutput_on
    [f,df]=feval(funame,x,varargin{:});
else
    f=feval(funame,x,varargin{:});
end
    
if t==-1
    f=f./(fm+l/2)*a;
    if GradOutput_on
        M_l = ones(length(x),1)*(fm+l/2);
        df=df./M_l*a;
    end
elseif t==0
else
    f=f(t);
    if GradOutput_on
        df=df(:,t);
    end
end