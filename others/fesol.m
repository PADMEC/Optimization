  function varargout=fesol(area,props,fext,glb,tipo)
%
% Realisa Análise via o MEF
% [u,sig,esf,vol,du,Lambd,D,freq]

    ang = props(1,:);
	comp = props(2,:);
	els = props(3,:); 
	ro = props(4,:);
if nargin<5
    tipo=[1 0 0];
end

global mat_curv
[sigi0,E0]=strain2stress(mat_curv,1e-16);
E=comp'*0+E0;
        
if tipo(1)>0 |tipo(2)>0
    
    %
    % Matriz de rigidez global:
    %
    kg = keglb(area,comp,E,ang,glb);
    %
    % Deslocamentos:
    %
    u = kg\fext;
    nnu=find(isnan(u));
    u(nnu)=0;
end
if tipo(1)>0

    %
    % Energia:
    en=fext'*u;

    %
    % Esforços:
    esf = tresf(area,comp,els,ang,glb,u);
    %
    % Tensões:
    sig = esf'./area';%sigma(area,comp,esf);
    %
    % Volume:
    vol = sum(area.*comp.*ro);%volume(area,comp);

    E=els;
    x=area;
    floc=sig./(-x'.*E'*pi/4./comp'.^2);
    varargout={u,sig,en,vol,floc};
end

%flamb=1;
du=[];D=[];
Lambd=zeros(1,length(fext));
freq=zeros(1,length(fext));

if tipo(2)>0
%floc=area.*els*pi/4./comp.^2;
%if nargout>7
%    [Lambd, du,uf,sigf,esff]=flambsol(esf,props,glb,kg,area,fext,u);
%    varargout={u,sig,en,vol,floc,du,Lambd,uf,sigf,esff};
    [Lambd, du]=flambsol(esf,props,glb,kg);
end
if tipo(3)>0

    [freq, D]=vibsol(props,glb,kg,area,ro);

end

varargout={u,sig,en,vol,floc,du,Lambd,D,freq};