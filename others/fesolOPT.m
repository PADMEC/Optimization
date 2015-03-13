  function varargout=fesolOPT(area,props,fext,glb,tipo,ROM,Zd)
  %
  % TEST LINEAR (u/du): RBM - SPACE
%
% Realisa Análise via o MEF
% [u,sig,esf,vol,du,Lambd,D,freq]

    ang = props(1,:);
	comp = props(2,:);
	els = props(3,:); 
	ro = props(4,:);

if nargin<7
    ROM=0;
end
if nargin<5
    tipo=[1 0 0];
end

global mat_curv
[sigi0,E0]=strain2stress(mat_curv,1e-16);
E=comp'*0+E0;
        load aaa
        fext=R;
if tipo(1)>0 |tipo(2)>0
    
        if (ROM==1)
            %if iter==1
                kg = keglb(area,comp,E,ang,glb);
                Kn=Zd'*kg*Zd;
                Kn1=pinv(Kn);
            %end
            Rn=Zd'*fext;
            alf = Kn1*Rn;
            u = Zd*alf;
            
            u2 = kg\fext;
            disp(norm(u-u2)/norm(u2));
        elseif (ROM==0)
            %
            % Matriz de rigidez global:
            %
            kg = keglb(area,comp,E,ang,glb);
            %
            % Deslocamentos:
            %
            u = kg\fext;
        end
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