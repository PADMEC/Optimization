function [fob,gob,fre,gre]=optsol(props,fext,glb,link,tobj,tres,xvalu,ndvab,lpdva,perturb)

        global apr  para  fs gkx itp
        global modos
        
% Fornece valores para as funcões objetivo e restrição e suas derivadas
%
%  atualiza todas as variáveis primárias
%
%	pega  nas variaveis   , "props"  alguns paramentros adicionais
%
    ang  = props(1,:);
	comp = props(2,:);
	els  = props(3,:);
%
fmof=0;
if tobj == 0
    global mo
    fmof=mo(1,3);
    if fs==3|fs==5
        %%%   NBI   %%%
        global t
        t=xvalu(end);
        xvalu(end)=[];
    end
end
    %area = getval(comp,xvalu,link);	
    area = xvalu(link(:,1)).*link(:,2)';	
%
%   Faz a solução via o MEF para avaliação das funções.
%
% 
%   Faz a análise de sensibilidades para avaliação dos gradientes das funções.
% 
%
%   guarda na variavel  ' 'par' alguns paramentros adicionais
% 			
    par = [ perturb ; ndvab];
%
% Armazena todas as funções e seus gradientes
%
    nmdfl=1;
    
    dlamb=zeros(ndvab,nmdfl);
    lamb=zeros(nmdfl,1);
    D=[];
    freq=[];
    [u,sig,en,vol,floc,V,lamb,D,freq]=fesol(area,props,fext,glb,modos);
    [den,du,dsig,dvol,dfloc,dlamb]=sensol(area,props,par,fext,glb,link,lpdva,sig,u,vol,floc,gkx,modos,V,lamb);

%if (sum(itp~=2)==tres)&(tobj~=3)&(fmof~=1)
%    [u,sig,en,vol,floc]=fesol(area,props,fext,glb);
%    
%    [du,dsig,dvol,dfloc]=sensol(area,props,par,fext,glb,link,lpdva,sig,u,vol,floc,gkx);
%else
%
%    %nl=ceil(length(lamb)/4);
%end
%disp(floc(1:10)')
    [fre,gre,fob,gob]=totim(tres,tobj,vol,dvol,en,den,sig,dsig,u,du,floc,dfloc,ndvab,area,props,lamb(1:nmdfl),dlamb(:,1:nmdfl));    

