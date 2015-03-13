function [f,df]=nbifun(x,funame,coname,tp,fm,l,a,n,par,varargin)
global GradOutput_on
df=[];
if tp==-1
    f=-x(end);
    df=x*0;df(end)=-1;
elseif tp==0
    %f=feval(funame,x,varargin{:});
    if GradOutput_on
        [f,df]=feval(funame,x,varargin{:});
        df(end+1,:)=f'*0;
    else
        [f]=feval(funame,x,varargin{:});
    end
else
    %f=feval(funame,x,varargin{:});
    if GradOutput_on
        [f,df]=feval(funame,x,varargin{:});
        f=f(tp);
        df=df(:,tp);
    else
        [f]=feval(funame,x,varargin{:});
        f=f(tp);
    end
end