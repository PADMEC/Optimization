function [Pt,n,ff,Pp,tp,Vr,xvr,P0]=projctpareto(ft,fs,minf,mnf)

    ff=ft;
    mnf=abs(mnf);
    mnf(mnf==0)=1;
    n_dim=length(fs);
    %for i_dim=1:n_dim
    f=fs{n_dim};
    [ndiv,nfun,npar]=size(f);
    
    for i=1:nfun
        ff(:,i)=(ff(:,i)-minf(i))/mnf(i);
    end
    P0=ff(end,:)';
    
%     Vr=[ff(2,:)'-ff(1,:)' ff(3,:)'-ff(1,:)'];
%     if norm(Vr(:,1))<1e-3
%         Vr=[f(end,:,1)'-f(end,:,2)' ff(3,:)'-ff(2,:)'];
%     end
%     if norm(Vr(:,2))<1e-3
%         Vr=[ff(2,:)'-ff(3,:)' f(end,:,1)'-f(end,:,3)'];
%     end
    for i=1:nfun-1
        Vr(:,i)=ff(end,:)-ff(i,:);
    end
    
    %NORMAL VECTOR GENERATE
    global mP
    Pt=null(Vr');
    if rank(Vr,1e-6)<(nfun-1)
        for ipar=1:npar
            for ifun=1:nfun
            mP(ipar,ifun)=mean((f(:,ifun,ipar)-minf(ifun))/mnf(ifun));
            end
        end
        for i=1:npar-1
            Vr(:,i)=mP(end,:)-mP(i,:);
        end

        Pt=null(Vr');
    end
        
    %Vr(:,2)=null([Vr(:,1) Pt]');
    %Vr(:,1)=Vr(:,1)/norm(Vr(:,1));
    
    n=-abs(Pt.*mnf')'/ndiv;
    
    %%%%%%%%%%%
    %PROJECTION
    %P0+Vr*abp+Pt*t=pf
    %[Vr Pt]*abt=pf-P0
    %A*abt=B
    pf=[];
    for i_dim=1:n_dim
        f=fs{i_dim};
        [ndiv,nfun,npar]=size(f);
        for ipar=1:npar
            pf=[pf; f(:,:,ipar)];
        end
    end
    for i=1:nfun
        pf(:,i,:)=(pf(:,i,:)-minf(i))/mnf(i);
    end
    [Pp,xvr,tp]=PointProj(Vr,P0,pf);
    
