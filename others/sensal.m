% sensibilidades = metodo direto semi-analitico

  function [ugra,sigra,volgra]=sensal(area,props,par,fext,glb,link,lpdva,sig,u,vol)
%####################################################################
%# Realisa Análise de sensibilidades usando o metodo direto semi-analitico	   #
%####################################################################
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
%
% Atualização das variáveis
%
    for jelem = 1:nelem
	    ivp = link(jelem,1);
		scal = link(jelem,2);
	    arean(jelem) = xpdva(ivp)*scal;
    end
%
% Calcule novos valores p/ as funcoes:
%
   [un,sign,esfn,voln] = fesol(arean,ang,comp,els,fext,glb);
%
% Cálculo dos gradientes via diferencas finitas
%
%	   ----volume----
%
%
    volgra(idvab) = (voln - vol)/perta; 
%
%	 ----deslocamentos---
% 
%   
    kgn=keglb(arean,comp,els,ang,glb);
    fextm(:,idvab)=-((kgn-kg)/perta)*u;
%
%	volta ao valor original
%
    xpdva(idvab) = area(ielem);
%
%
%
end
%fextm;
ugra=kg\fextm;

%
%	 ----tensões---
%  
for idvab = 1:ndvab
   ielem = lpdva(idvab);
      
   for jelem = 1:nelem
      kdvab = link(jelem,1); 
      kelem = lpdva(kdvab);
            
     % if kelem == ielem
         Kelem = [-cos(ang(jelem)) -sin(ang(jelem)) cos(ang(jelem)) sin(ang(jelem))];
         id = find(glb(:,jelem));
         gb = glb(id,jelem);
         duelem = ugra(gb,idvab);
         den = Kelem(id)*duelem;
         dkl = els(jelem)/comp(jelem);
         sigra(idvab,jelem) = dkl*den;
    %  else  
    %     sigra(idvab,jelem) = 0;
    %  end;
      
   end;
end;
ugra=ugra';
