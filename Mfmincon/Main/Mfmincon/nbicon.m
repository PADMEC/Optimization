function [R,C,dR,D]=nbicon(x,funame,coname,tp,fm,l,a,n,par,varargin)
global GradOutput_on
%dR=[];
C=[];
D=[];
if tp==-1
    t=x(end);
    x(end)=[];

    if GradOutput_on
        [R,~,dR]=feval(coname,x,varargin{:});
        [f,df]=feval(funame,x,varargin{:});
        dR(end+1,:)=R'*0;
        df(end+1,:)=f'*0;
        df=df(:,par);
        dt=x*0;dt(end+1)=1;
        dres=dt'*n'-df;
        dR=[dR -dres];
        %dR=[dR dres];
    else
        [R]=feval(coname,x,varargin{:});
        [f]=feval(funame,x,varargin{:});
        dR=[];
    end
    
    f=f(par)-fm(:)';
    res=a+t*n-f';
    R=[R(:); -res];
   % R=[R; res];
else
    if GradOutput_on
        [R,~,dR]=feval(coname,x,varargin{:});
    else
        [R]=feval(coname,x,varargin{:});
        dR=[];
    end
end