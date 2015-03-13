function [f,df]=mmxfun(x,funame,coname,t,mf,l,a,n,par,varargin)


global GradOutput_on
if GradOutput_on
    [f,df]=feval(funame,x,varargin{:});
else
    f=feval(funame,x,varargin{:});
end

if t==-1
    f=(f-mf)./l;
    [f,imax]=max(f.*a');
    if GradOutput_on
        M_l = ones(length(x),1)*l;
        df=df./M_l;
        df=df(:,imax)*a(imax);
    end
    
elseif t>0
    f=f(t);
    if GradOutput_on
        df=df(:,t);
    end
    
end