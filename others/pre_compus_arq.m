function [arq,area,props,fext,glb,tipo,nel,mdf,mdv]=pre_compus_arq(arq)

%
% Entrada de Dados relacionado com a Geometria dos elementos
%
% forneça as coordenadas dos nos
%load  C:\MATLAB7\work\arq0

global X Y Z conect x0 link Truss3D


if isfield(arq,'noz')
    % Spartial 3D Truss
    Truss3D=1;
    Z=arq.noz';
    glpn=3;
else
    % 2D Truss
    Truss3D=0;
    glpn=2;
end

X=arq.nox';Y=arq.noy';conect=arq.barra;
nel=size(arq.barra,1);nno=length(arq.nox);

% Calculando o comprimento e os angulos das barras
if Truss3D
	[cosn,comp] = compang3D(X,Y,Z,conect) ;
else
    [cosn,comp] = compang  (X,Y,conect) ;
end

% forneça a area da seçao tansversal de cada elemento
%if isfield(arq,'fp')
%    if isfield(arq,'x')
%        if size(arq.x,2)==nel
%            area=arq.x;
%        else
%            area=arq.x(arq.arex);
%        end
%    else
%        area = arq.x0*ones([1,nel]);
%    end
%else
%    area = arq.x0*ones([1,nel]);
%end

if size(arq.arex,1)==1&&size(arq.arex,2)>1
    arq.arex=arq.arex';
end
    
    
    

if isfield(arq,'x0')
    
    if ~isfield(arq,'arex')
        nelm=size(arq.barra,1);
        arq.arex=1:nelm;
    end
    if size(arq.x0,2)~=max(arq.arex(:,1))
        arq.x0=arq.x0(1)*ones(1,max(arq.arex(:,1)));
    end
else
    arq.x0=1;
    arq.x0=arq.x0*ones(1,max(arq.arex(:,1)));
end

if size(arq.x0,2)==max(arq.arex(:,1))
    %x0=arq.x0;
    area = arq.x0(arq.arex(:,1));
elseif length(arq.x0)==length(arq.arex)
    area = arq.x0;
elseif length(arq.x0)==1
    x0=arq.x0*ones(1,max(arq.arex(:,1)));
    area = x0(arq.arex(:,1)); 
else
    area = ones(1,max(max(arq.arex)));
end

if length(arq.E)~=nel
    els=arq.E*ones(1,nel);
else
    els = arq.E;
end


% forneça o peso especifico do material
if length(arq.ro)~=nel
    ro=arq.ro*ones(1,nel);
else
    ro = arq.ro;
end

% guarda as propriedades do material

props = {cosn;comp;els;ro};

%Dados relacionado com Graus de liberdade globais da 
%estrutura construindo a matriz ID
%			0 => Grau de liberdade livre (desconhecido)
%			1 => grau de liberdade restrito (conhecido)
nax=size(arq.ax,2);
nay=size(arq.ay,2);

    
ID = zeros(nno,glpn);
ID(arq.ax,1)=ones(nax,1);
ID(arq.ay,2)=ones(nay,1);
if Truss3D
    naz=length(arq.az);
    ID(arq.az,3)=ones(naz,1);
end

% Degrees of Freedom
global ngln
ngln=1:nno*glpn;
igl_out=[(arq.ax-1)*glpn+1 (arq.ay-1)*glpn+2];
if Truss3D
    igl_out = [igl_out (arq.az-1)*glpn+3];
end
ngln(igl_out)=[];


%Construindo a matriz dos graus de liberdade

	glb=constroi(ID,conect);

%
% Entrada de Dados relacionado com as Cargas externas
%
fext=zeros(glpn*nno,1);% -nax-nay
fext(1:glpn:end)=arq.px;
fext(2:glpn:end)=arq.py;
if Truss3D
    fext(3:glpn:end)=arq.pz;
end
fext(igl_out)=[];
    
%		1---MEF ANALISE
%
%%%%%%%%%%%%%%%%%%%%% P/ O COPILADOR (.exe) %%%%%%%%%%%%%%%%%%%%%
aaa=1;
if aaa==0
    mycon(x,xol,fext,glb,link);
    myfun(x,xol,fext,glb,link);
    myconrm2(x,xol,fn,N,Z);
    myfrm2(x,xol,fn,N,Z);
    axesCratFcn(hObject,tp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mdf=arq.modos(2);
mdv=arq.modos(3);
tipo=arq.modos;
%arq.mdv=1;
%tipo=[arq.mdf arq.mdv];
disp('Areas')
disp(arq.x0);

x0=arq.x0;
link=[arq.arex(:),arq.arex(:)*0+1];
end