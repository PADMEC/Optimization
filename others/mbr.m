function arq=mbr(arq,plt,SN)

% Programa para analise de treliças planas pelo metodo
% dos deslocamentos
%
% opção desejada:   1==análise via o MEF,
%				        2==análise de sensibilidades
%                   3==otimização.
%
    global para
    global  xol
    global dvol
    global lpdva
    global props
    global modos
    global link glb els ang ro comp
%
  %createtrussge
  tic
 arq2rbm(arq);
 
 Precomputime=toc
  

%
%reduced basis version
%input data for the basis parameters
%

  
%SN=gridsamp([.5 .5 .5;8 7 7],2);
N= size(SN,1);
nel=size(arq.barra,1);
em=[];cm=[];d=[];
try
if length(arq.em)~=nel
    em=arq.em*ones(1,nel);

else
    em=arq.em;
end
arq.cm=-abs(arq.cm);
if length(arq.cm)~=nel
    cm=arq.cm*ones(1,nel);
else
    cm=arq.cm;
end
end

if size(arq.x0,2)~=max(arq.arex)
    x0=arq.x0*ones(1,max(arq.arex));
else
    x0=arq.x0;
end
munew = x0;   
  

  load trussdata.mat


% Execute a opção desejada:
%		3--- OTIMIZAÇ+O
%

             ndvab = size(SN,2);
             xpdva = munew;
             area=munew(link(:,1));
             
 lpdva=ones(1,ndvab);
 for i=1:ndvab
    iv=find(link(:,1)==i);
    lpdva(i)=iv(1);
 end
  
  ngl=size(K1{1},1);
  nel=size(link,1);
  
  [fde,fte,tpres]=fobrest(arq,ngl,area,els,ro,comp,cm,em,d);

        %tobj=2;%    1---VOLUME  2----ENERGIA     3---Dis.    4--Disl. Especifc.
        %tpres=1;%    1---PESO    2----TENSAO     3---Tens e Desloc   4---Desloc
        vlb =  arq.xmin * ones(1,ndvab);
        vub =  arq.xmax* ones(1,ndvab);

        itera = 0;
		 para = [1e-6;ndvab;arq.tobj;tpres;itera];
	     %props = [ang;comp;els;ro];
         modos=arq.modos;

         % Criaçao de Z(desl. das amostras) e valores iniciais
tic

    [fn,kni,Z,kvi,Z0,kdi,kmni,Zd]  =  feoffline2(SN,modos);%No. de modos de flambagens

Off_Line=toc

    [alph,arq.idis,arq.ites,arq.ien,arq.ivol,floc,alfv,du,Lambd,k0n] = feonline2(munew,modos,fn,kni,Z,kvi,Z0,kdi,kmni,Zd);

    if arq.tpres(2)==1|(arq.tobj==0&arq.mo(3)==1)

    else
        Lambd=-1;
        kvi=[];
        Z0=[];
    end


arq.ilmbd=Lambd(1);

       %[dsn]=sensolrm2(alph,munew,ndvab,N,fn,kn,kni);
       
       
                        %fs=3;
if arq.tobj==0      %  MULTI-OBIJETIVO
    
    %arq=morbm(arq,plt,munew,vlb,vub,xpdva,'MObjrbm',fn,kni,Z,kvi,Z0);
    arq=moptruss(arq,plt,munew,vlb,vub,xpdva,'MObjrbm',2,fn,kni,Z,kvi,Z0);
    %arq=moptrussrbm(arq,plt,x0,vlb,vub,xpdva,fext,glb,link);
else 
       

    %       
    %  Uses optimization routine (fmincon)

    %                           
         options=optimset('LargeScale','off');
           options=optimset(options,'GradObj','on','GradConstr','on','Display','off');
           options=optimset(options,'Display','iter','DerivativeCheck','off');%,'MaxIter',12);%,'MaxTime',9);

           n=length(xpdva);
    %       
    %             final solution
    %
    itp=find(arq.tpres);
    tic;
    [arq.x,arq.ffobj,arq.Converg,arq.Resumo]=fmincon('myfrm2',munew,[],[],[],[],vlb,vub,'mycrm2',options,xpdva,fn,kni,Z,kvi,Z0);
    arq.t=toc; 
    disp('BASES :')
    disp(SN)
    arq.Resumo
    disp('X =')
    disp(arq.x)
    arq.x=arq.x(arq.arex);
    disp('Função Objeto')
    disp(arq.ffobj)


    %  
    %  [arq.fdis,arq.ftes,esf,arq.fvol,alph,fn,kn,kni,Z]=fesolrm2(arq.x,SN,N);

end
    props = [ang;comp;els;ro];
    [arq.fdis,arq.ftes,arq.fen,arq.fvol,floc,du,Lambd,D,freq]=fesol(arq.x,props,F,glb,[1 0 0]);

    disp('Lambda = ');disp(Lambd(1))%:ceil(nglb/2))'
    arq.flmbd=Lambd(1);

      disp('Exato')
      switch arq.tobj
          case 1
              disp(arq.fvol)
          case 2
              disp(arq.fen)
          case 4
              disp(norm(arq.fdis))
          case 5
              disp(max(abs(arq.fdis)))
          case 6
              disp(max(abs(arq.ftes)))
          case 7
             % disp(arq.fdis(arq.idm))
      end

    disp('tempo de otimizacao')
    disp(arq.t)
