function arq2rbm(arq)
% Programa para analise de treliças planas pelo metodo dos deslocamentos
%
% opção desejada:   1==análise via o MEF,
%				    2==análise de sensibilidades
%                   3==otimização.
%
%

%
% Entrada de Dados relacionado com a Geometria dos elementos
%
% forneça as coordenadas dos nos

global link glb els ang ro comp

%[barra,ax,ay,nox,noy,px,py,arex,E,em,cm,dd,ro,x00,xmin,xmax]=feval(fun);
X=arq.nox';Y=arq.noy';conect=arq.barra;
nel=size(arq.barra,1);nno=size(arq.nox,2);
if size(arq.x0,2)~=max(arq.arex)
    x0=arq.x0*ones(1,max(arq.arex));
else
    x0=arq.x0;
end
try
if length(arq.em)==1
    arq.em=arq.em*ones(1,nel);
    arq.cm=arq.cm*ones(1,nel);
end
end
%arq.d=3;
% Calculando o comprimento e os angulos das arq.barras
	[ang,comp] = compang(X,Y,conect) ;

% Forneça o modulo de elasticidade de cada arq.barra
%E = 10^4;
%arq.E =.207e9;
%arq.E = 30000;
els = arq.E*ones([1,nel]);

% forneça o peso especifico do material

ro = arq.ro*ones([1,nel]);

% guarda as propriedades do material

props = [ang;comp;els;ro];

nax=size(arq.ax,2);
nay=size(arq.ay,2);
ID = zeros([nno,2]);
ID(arq.ax,1)=ones(nax,1);
ID(arq.ay,2)=ones(nay,1);
glb=constroi(ID,conect);

fext=zeros(2*nno,1);% -nax-nay
fext(1:2:end)=arq.px;
fext(2:2:end)=arq.py;
fext([arq.ax*2-1 arq.ay*2])=[];

F=fext;
link=[arq.arex(:) ones(nel,1)];
ndv=max(arq.arex);
nb=nel;
n=nno;

for i=1:nb
    dx=arq.nox(arq.barra(i,1))-arq.nox(arq.barra(i,2));
    dy=arq.noy(arq.barra(i,1))-arq.noy(arq.barra(i,2));
    l(i)=(dx^2+dy^2)^.5;
    c(i)=dx/l(i);
    s(i)=dy/l(i);    
    imr(:,:,i)=[c(i)^2 c(i)*s(i);c(i)*s(i) s(i)^2];
end
for i=1:ndv
    bdv=find(arq.arex==i);
    dvol(i,1)=sum(l(bdv).*ro(bdv));
end

k=arq.E./l;
%%%%%%%  Matriz Rigidez   %%%%%%%%%%%
mr=zeros(n*2,n*2,ndv);

%   Diagonal
for i=1:n;
    [b d]=find(arq.barra==i);
    k1=i*2-1;
    for k2=1:size(b)
    mr(k1:k1+1,k1:k1+1,arq.arex(b(k2)))=mr(k1:k1+1,k1:k1+1,arq.arex(b(k2)))...
        +k(b(k2)).*imr(:,:,b(k2));
    end
end

% Matriz Rigidez
for i=1:nb
    k1=arq.barra(i,1)*2-1;
    k2=arq.barra(i,2)*2-1;
    mr(k1:k1+1,k2:k2+1,arq.arex(i))=-k(i).*imr(:,:,i);
    mr(k2:k2+1,k1:k1+1,arq.arex(i))=-k(i).*imr(:,:,i);
end

%%%%%%%%%   Matriz Liberd    %%%%%%%%

A=[arq.ax*2-1 arq.ay*2];
na=size(A,2);
j=0;
for i=1:2*n
    d=find(A==i);
    if isempty(d)
        j=j+1;
        Ad(j)=i;
        %Pd(j)=P(i);
    end
end

mr=mr(Ad,Ad,:);

for i=1:ndv
    K1{i}=sparse(mr(:,:,i));
    arei=find(arq.arex==i);
    
    %Matriz de Massa unit
    Km1{i} = mmglb(arq.arex==i,comp,glb,ro);

end

save trussdata.mat K1 F comp link els glb ang dvol ro Km1
% xmin=arq.xmin;
% xmax=arq.xmax;
% 
% save ENTRADAS-1800LTS X Y ID glb link fext xmin xmax els props x0
