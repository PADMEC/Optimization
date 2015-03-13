function [R,C,dR,D]=mmxcon(x,funame,coname,t,mf,l,a,n,par,varargin)

C=[];D=[];
global GradOutput_on
if GradOutput_on
    [R,~,dR]=feval(coname,x,varargin{:});
else
    R=feval(coname,x,varargin{:});
end

%R=R(:)';
if t==-1
    if GradOutput_on
        [f,df]=feval(funame,x,varargin{:});
    else
        f=feval(funame,x,varargin{:});
    end
    
    f=(f-mf)./l;
    fa=f.*a';
    [fmax,km]=max(fa); 
    %nf=length(f);
    %imin=1:nf;imin(km)=[];
    fa(km)=[];
    R=[R fa-fmax];
    
    if GradOutput_on
        M_l = ones(length(x),1)*(l./a');
        dfn=df./M_l;
        
        nf=length(f);
        imin=1:nf;imin(km)=[];
        dR=[dR (dfn(:,imin)-(dfn(:,km)*ones(1,nf-1)))];
    end
end