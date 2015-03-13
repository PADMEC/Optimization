% sensibilidades = metodo analítico

  function [den,ugra,sigra,volgra,gfloc,lambgrad]=senal(area,props,par,fext,glb,link,lpdva,sig,u,gkx,modos,V,Lambd)

%####################################################################
%# Realisa Análise de sensibilidades usando o metodo analítico      #
%####################################################################
%Obs:volgrad é o gradiente do PESO
%
%  Entrada de dados 
%
%
%	pega  nas variaveis par  , "props"  alguns paramentros adicionais
%
%
    ang = props(1,:);
    comp = props(2,:);
	els = props(3,:); 
	ro = props(4,:);
%
    perturb =  par(1);
	ndvab = par(2);
%	
%% dimensoes:
%
[m,nelem]=size(comp);
 nglb = max(max(glb));
%
% Inicializa os gradientes:
%
ugra = zeros(ndvab,nglb);
sigra = zeros(ndvab,nelem);
volgra = zeros(ndvab,1);
xpdva =   zeros(1,ndvab);
%
% vetor das variáveis primárias
%
 for idvab = 1:ndvab
     iel = lpdva(idvab);
     xpdva(idvab) = area(iel);
 end 

kg=keglb(area,comp,els,ang,glb);
esf=sig'.*area;
KG=kgeom(esf, comp, ang, glb);

E=els;
x=area;

c=cos(ang);
s=sin(ang);
Ec=E./comp;

% 
% Loop sobre todas as variaveis primarias:
%
for idvab = 1:ndvab
%
% Perturbe isoladamente cada variavel de projeto:
%
    ielem = lpdva(idvab);
    xpdva(idvab) = area(ielem)*(1 + perturb);
    perta = area(ielem)*perturb;
    
    idv=find(link(:,1)==idvab);
%
%	 ----deslocamentos---
% 
% 
    fextm=-gkx(:,:,idvab)*u;
    ugra(idvab,:)=kg\fextm;
    
    den(idvab) = ugra(idvab,:)*fext;
    

%
% Atualização das variáveis
% 
     du=ugra(idvab,:);
     %Deslocamento=0=>glb=0
     glb(find(glb==0))=size(du,2)+1;
     du(end+1)=0;
%     
     sigra(idvab,:)=c.*(du(glb(3,:))-du(glb(1,:))).*Ec + s.*(du(glb(4,:))-du(glb(2,:))).*Ec;
% dsig(idvab,:)=1;
%     
     glb(find(glb==size(du,2)))=0;
     du(:,end)=[];

    for jelem = 1:nelem
	    ivp = link(jelem,1);
	    scal = link(jelem,2);
		arean(jelem) = xpdva(ivp)*scal;
        
%
%	 ----tensões---
%      
            
        %if kelem == ielem
%         ae = [-cos(ang(jelem)) -sin(ang(jelem)) cos(ang(jelem)) sin(ang(jelem))]; 
%         id = find(glb(:,jelem));
%         gb = glb(id,jelem);
%         duo = ugra(idvab,gb)';
%         den = ae(id)*duo;
%         dkl = els(jelem)/comp(jelem);
%         sigra(idvab,jelem) = dkl*den;
    end
%
%	 ----tensões criticas---
%      
        gfloc(idvab,:) = sigra(idvab,:)./(-x.*E*pi/4./comp.^2);
        gfloc(idvab,idv)= gfloc(idvab,idv)+sig(idv)'./(E(idv)*pi/4./comp(idv).^2)./x(idv).^2;
    
%
%	   ----volume----
%
    peselem=comp.*ro;
    volgra(idvab)=sum(peselem(find(link(:,1)==idvab)));
    
    lambgrad=zeros(ndvab,2);
    if modos(2)>0
%
%        ----Flambagem Global---
%
%        ----Esforço da Vab---
%
        desfi=sigra(idvab,:).*area;
        desfi(idv)=desfi(idv)+sig(idv)';

        gk0=kgeom(desfi, comp, ang, glb);
        for i=1:length(Lambd)
            normlz=abs(V(:,i)'*KG*V(:,i));
            lambgrad(idvab,i)=V(:,i)'*(gkx(:,:,idvab)+Lambd(i)*gk0)*V(:,i)/normlz;
        end
    end
%
%	volta ao valor original
%
    xpdva(idvab) = area(ielem);
%

end

