function [dsn,dsig,du,gfloc,lambgrad]=seonline(munew,kn,props,glb,link,modos,sig,alpha,fn,kni,Z,kvi,Z0,k0n,alfv,Lambd);
%
%grad(F,X)=[dF3/dX1   dF2/dX1  dF3/dX1  dFnf/dX1]
%          [dF3/dX2   dF2/dX2  dF3/dX2  dFnf/dX2]
%          [dF3/dX3   dF2/dX3  dF3/dX3  dFnf/dX3]
%          [dF3/dXnx dF2/dXnx dF3/dXnx dFnf/dXnx]
%   
%  online stage for sensitivity analysis:
%
%

 % 1) pseudo load vecto fn*:
 %   fn* = dfn/dxi - dkn/dxi*alpha
 %
 %   dkn/dxi:
%
%loop over dv  
%

area=munew(link(:,1));

%props = [ang;comp;E;ro];
ang = props(1,:);
comp = props(2,:);
E = props(3,:); 
ro = props(4,:);

ndvab=size(kni,3);
if 1==2%abs(det(kn))>1e-3
    invKn=kn^-1;
    for i=1:ndvab
        f1=-kni(:,:,i)*alpha;
        dalpha=invKn*f1;
        du(i,:)=Z*dalpha;
        dsn(i) = dalpha'*fn;
    end
else
    for i=1:ndvab
        f1=-kni(:,:,i)*alpha;
        dalpha=kn\f1;
        du(i,:)=Z*dalpha;
        %du(i,:)=Z(:,1)';
        dsn(i) = dalpha'*fn;
    end
end

x=area;
c=cos(ang);
s=sin(ang);
ndb=length(c);
Ec=E./comp;
% 
% Loop sobre todas as variaveis primarias:
%
for j=1:ndvab

    idv=find(link(:,1)==j);    
    %
    %	 ----tensões---
    %
    glb(find(glb==0))=size(du,2)+1;%Deslocamento=0=>glb=0
    du(:,end+1)=0;
    %tic
    
    dsig(j,:)=c.*(du(j,glb(3,:))-du(j,glb(1,:))).*Ec + s.*(du(j,glb(4,:))-du(j,glb(2,:))).*Ec;
%dsig(j,:)=1;
    %for i=1:ndb
    %    dui=du*E(i)/comp(i);
    %    dsig(j,i)=c(i)*(dui(j,glb(3,i))-dui(j,glb(1,i)))...
    %        +s(i)*(dui(j,glb(4,i))-dui(j,glb(2,i)));
    %end
    %toc
    glb(find(glb==size(du,2)))=0;
    du(:,end)=[];
    
    %dsig(j,:)=dsig(j,:);
    
    %	 ----tensões criticas---
    gfloc(j,:) = dsig(j,:)./(-x.*E*pi/4./comp.^2);
%gfloc(j,:)=1;
    gfloc(j,idv)= gfloc(j,idv)+sig(idv)./(E(idv)*pi/4./comp(idv).^2)./x(idv).^2;

    lambgrad=zeros(ndvab,size(Z,1));
    
    if modos(2)>0
    
%        ----Esforço da Vab---
%
        desfi=dsig(j,:).*area;
        %desfi=desfi-desfi;
        desfi(idv)=dsig(j,idv).*area(idv)+sig(idv);
%
%        ----Flambagem Global---
%   
        gk0=kgeom(desfi, comp, ang, glb);
        for imdf=1:length(Lambd)
            gk0n = Z0(:,:,imdf)'*gk0*Z0(:,:,imdf);
            normlz=alfv(:,imdf)'*k0n(:,:,imdf)*alfv(:,imdf);
            lambgrad(j,imdf)=alfv(:,imdf)'*(kvi{j}(:,:,imdf)+Lambd(imdf)*gk0n)*alfv(:,imdf)/normlz;
            %%Conferencia
            %load trussdata.mat
            %V=Z0(:,:,imdf)*alfv(:,imdf);
            %KG=kgeom(sig.*area, comp, ang, glb);
            %normlz=V'*KG*V;
            %lambgrad(j,imdf)=V'*(K1{j}+Lambd(imdf)*gk0)*V/normlz;

            %kvi=>{n da variavel=j}(gral de liberdade(N),n da amostra(N), n do modo d falmbagem)
        end
    end
    if modos(3)>0
        
        disp('fazer gradiente de vibrações no RBM (seonline)')
        %OPS!
    end
end