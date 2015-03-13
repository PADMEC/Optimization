function     [f,df]=nncfun(x,funame,coname,t,fm,l,xn,N,par,varargin)

global GradOutput_on
if GradOutput_on
    [f,df]=feval(funame,x,varargin{:});

    if t==-1
        M_l = ones(length(x),1)*l;
        df=df(:,par)./M_l;
        df=df(:,end);
    elseif t>0
        df=df(:,t);
    end
else
    f=feval(funame,x,varargin{:});
end

if t==-1
    f=(f(par)-fm)./l;
    f=f(end);
elseif t>0
    f=f(t);
end