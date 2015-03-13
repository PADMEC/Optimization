function varargout= feonline2(munew,tipo,fn,kni,Z,kvi,Z0,kdi,kmni,Zd)
%
%
%  online stage for FE analysis:
%
%
% first build KN matrix per turn per each new parameter mu
%
% 
   %load trussdata.mat
global link glb ang comp els ro 
%LargTruss
%munew(3)=1
 area=munew(link(:,1));
%
% volume
%
vol =sum(comp.*area.*ro);

if tipo(1)>0
    
    kn=kni(:,:,1)*0;
    for i=1:size(kni,3)
        kn = kn+ munew(i)*kni(:,:,i);
    end

%kn=sparse(kn);
   %kn = munew(1)*kni;
%
% then solve Kn*alph = Fn
%
   alph =  kn\fn;

%
% finally calculates outputs
%
  sn = alph'*fn;

    u=Z*alph;
    %u=Z(:,1);

    esf = tresf(area,comp,els,ang,glb,u);
%
% Tensões:
%
    sig = esf./area;
    
 E=els;
 x=area;
 floc=sig./(-x.*E*pi/4./comp.^2);
%floc=1;

end
alfv=[];du=[];k0n=[];D=[];
Lambd=zeros(1,size(Z,1));
freq=zeros(1,size(Z,1));

%FLAMBAGEM
if tipo(2)>0

    du=[];
    Lambd=[];
    alfv=[];
    %nmdf=2;%SO O PRIMEIRO E O SEGUNDO MODOS DE FLAMBAGEM
    nmdf=size(Z0,3);
    k0 = kgeom(esf,comp,ang,glb);
    for imdf=1:nmdf
        k0n(:,:,imdf)=Z0(:,:,imdf)'*k0*Z0(:,:,imdf);
        kv=kvi{1}(:,:,1)*0;
        for idv=1:size(kvi,2)
            kv = kv+ munew(idv)*kvi{idv}(:,:,imdf);
        end
        [av,D]=eig(kv,k0n(:,:,imdf));
        D=diag(D);
        %Modos de Flambagem
        [D,mdf]=sort(abs(real(D)));
        av=av(:,mdf);
        for iver=1:size(av,2);
            if norm((kv-D(iver)*k0n(:,:,imdf))*av(:,iver))/norm(kv*av(:,iver))<1e-4
                break
            end
        end
        
        alfv=[alfv av(:,iver)];
        V=Z0(:,:,imdf)*av(:,iver);
        du=[du full(real(V))];
        Lambd=[Lambd D(iver)];
    end

    %Confere
    %K=K1{1}-K1{1};
    %for i=1:size(K1,2)
    %    K = K+munew(i)*K1{i};
    %end
    %tic
    %[V,D]=eig(full(K),full(k0));
    %D=diag(D);
    %[Lambdr,mdf]=sort(abs(real(D)));
    %du1=full(real(V(:,mdf(1:nmdf))));
    %Lambd1=Lambdr(1:nmdf);
    %if abs((Lambd(2)-Lambd1(2))/Lambd(2))>0.1
    %    [Lambd1(1:2) Lambd']
    %end
    %toc
    
end

%VIBRAÇÃO
if tipo(3)>0

    du=[];
    freq=[];
    alfd=[];
    %nmdf=2;%SO O PRIMEIRO E O SEGUNDO MODOS DE FLAMBAGEM
    nmdv=size(Zd,3);
    for imdv=1:nmdv
    kdvi=zeros(size(kdi{1}));
    kmnvi=zeros(size(kmni{1}));
        %Matriz de rigidez projetada (Kdi) no modo de vib (imdv)
        for idv=1:size(kdi,2)
            %kdmv = kd(:,:,imdv) + munew(idv)*kdi{idv}(:,:,imdv);
            kdvi = kdvi + munew(idv)*kdi{idv}(:,:,imdv);
        end
        %kd(:,:,imdv)=kdmv;
        
        %Matriz de massa projetada (Kmni) no modo de vib (imdv)
        for idv=1:size(kmni,2)
            %kmnmv = kmn(:,:,imdv) + munew(idv)*kmni{idv}(:,:,imdv);
            kmnvi = kmnvi + munew(idv)*kmni{idv}(:,:,imdv);
        end
        %kmn(:,:,imdv)=kmnmv;
        
        %kg*du+Lambd*mg*du==Indet
        %Z'*kg*Z*ad+Lambd*Z'*mg*Z*ad==Indet
        %(kdmv+Lambd*kmnmv)*ad==Indet
        [ad,D]=eig(kdvi,kmnvi);
        L=diag(D);
        %Ordenando Frequencias de vib
        [L,mdv]=sort(abs(real(L)));
        ad=ad(:,mdv);%ordenando
        for iver=1:length(ad);
            adi=ad(:,iver);
            if norm((kdvi-L(iver)*kmnvi)*adi)/norm(kdvi*adi)<1e-3
                break
            end
        end
        
        alfd=[alfd ad(:,iver)];
        D=Zd(:,:,imdf)*ad(:,iver);
        du=[du full(real(D))];%Modos de Vibração
        freq=[freq (L(iver))^.5];
    end

    %Confere
    %K=K1{1}-K1{1};
    %for i=1:size(K1,2)
    %    K = K+munew(i)*K1{i};
    %end
    %tic
    %[V,D]=eig(full(K),full(k0));
    %D=diag(D);
    %[Lambdr,mdf]=sort(abs(real(D)));
    %du1=full(real(V(:,mdf(1:nmdf))));
    %Lambd1=Lambdr(1:nmdf);
    %if abs((Lambd(2)-Lambd1(2))/Lambd(2))>0.1
    %    [Lambd1(1:2) Lambd']
    %end
    %toc
    
end
    varargout={alph,u,sig,sn,vol,floc,kn,alfv,du,Lambd,k0n,D,freq};
