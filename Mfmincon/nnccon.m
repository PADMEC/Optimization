function [R,C,dR,D]=nnccon(x,funame,coname,t,fm,l,xn,N,par,varargin)

C=[];D=[];

global GradOutput_on F_iters Ps

if GradOutput_on
    [R,~,dR]=feval(coname,x,varargin{:});
else
    R=feval(coname,x,varargin{:});
end

if t==-1
    %fn(1)=(x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2-fm(1))/l(1);
    %fn(2)=(3*x(1)+2*x(2)-x(3)/3+0.01*(x(4)-x(5))^3-fm(2))/l(2);
    %fn(3)=(((x(3)+x(4))^2+1/(100+x(2)))-fm(3))/l(3);
    
    if GradOutput_on
        [f,df]=feval(funame,x,varargin{:});
        
        M_l = ones(length(x),1)*l;
        dfn=df(:,par)./M_l;
        dR=[dR dfn*N];
        
    else
        f=feval(funame,x,varargin{:});
    end
    
    f=f(par);
    fn=(f-fm)./l;
    R=[R (fn-xn')*N];
    
    if sum(par)==sum(1:length(f))
        %1 2 3...nfun?
        Ps(end+1,:)=xn'.*l+fm;
        F_iters(end+1,:)=f;
    end
end