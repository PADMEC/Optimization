function [f,df]=nbifun(x,funame,coname,tp,fm,l,a,n,par,varargin)
global GradOutput_on
df=[];
test_proj_modf=0;
if tp==-1
    ft=-x(end);
    df=x*0;df(end)=-1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test NBI distance
    if test_proj_modf
        x(end)=[];
        if GradOutput_on
            [f,df]=feval(funame,x,varargin{:});
            df(end+1,:)=f'*0;
        else
            [f]=feval(funame,x,varargin{:});
        end
        fn=f(par)-fm(:)';
        %F_project_into_n=
        %Ft=P0+[n'*(P-P0)/norm(n)]*n/normn
        %Ft=P0+[n'*(P-P0)/(n'n)]*n
        %Ft=P0+t*n
        %Dist_fn=norm(Ft-fn)


        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PROJECTING F INTO N VECTOR
        %Ft=P0+[nu'*(P-P0)]*nu
        %Ft=P0+[nu'*nu](P-P0)
%       t=n'*(fn'-a)/(n'*n);
%       Ft=a+t*n;

        %Normalized projection (l)
        Dist=(Ft'-fn)./l;
        lu=l'/norm(l);
        nnu=n/norm(n)./lu;
        tn=nnu'*((fn'-a)./lu)/(nnu'*nnu);
        Ftn=a+tn*nnu.*lu;
        Dist=(Ftn-fn')./l';
        ft=ft + Dist'*Dist;%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f=ft;% + Dist'*Dist;%
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