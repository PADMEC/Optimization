function arq=mbranalestrut(arq,SN,modflam,sens)

% Programa para analise de treliças planas pelo metodo
% dos deslocamentos
%
% opção desejada:   1==análise via o MEF,
%				        2==análise de sensibilidades
%                   3==otimização.
%
    global para
    global con
    global  xol
    global lpdva
    global props
%
  %createtrussge
 arq2rbm(arq);
  

%
%reduced basis version
%input data for the basis parameters
%

  
%SN=gridsamp([.5 .5 .5;8 7 7],2);
N= size(SN,1);
nel=size(arq.barra,1);
em=[];cm=[];d=[];

if size(arq.x0,2)~=max(arq.arex)
    x0=arq.x0*ones(1,max(arq.arex));
else
    x0=arq.x0;
end
  

%disp('BASES :')
%SN

area = x0(arq.arex);
arq.x=area;
%disp('BASES :')
%disp(SN)
%disp('X =')
%disp(x0)
arq.mdf=modflam;

if sens==1
    [t,ugra,sigra,volgra,den,flocgrad,lambgrad]=sensolrm2(x0,SN,arq.modos);

    if arq.modos(2)==1
        disp('Fator de Flambagem')
        disp(lambgrad)
    end
    disp('Deslocamento')
    %disp(ugra)
    disp('Tensao')
    %disp(sigra)
    disp('Volume')
    disp(volgra')
    disp('Energia')
    disp(den)
    disp('Tempo para Analise')
    disp(t)
else
    [fn,kni,Z,alph,arq.t,arq.idis,arq.ites,arq.ien,arq.ivol,floc,alfv,du,Lambd,D,freq]=fesolrm2(x0,SN,arq.modos);
    arq.fdis=arq.idis;
    arq.ftes=arq.ites;
    arq.fen=arq.ien;
    arq.fvol=arq.ivol;

        disp ('Compliance')
        disp (arq.fen)
        
    %FLAMBAGEM
    if arq.modos(2)==0
        
        lbd=1;
        arq.flmbd=0;
        arq.ilmbd=0;
        arq.escflamb=0;
    else
        disp ('Lambd')
        Lambd=Lambd(arq.modos(2));
        arq.escflamb=1;
        arq.D=D(:,arq.modos(2));
        disp(Lambd)
        %arq.x=area;%/Lambd;
        %load trussdata.mat
        %esff = tresf(arq.x,comp,els,ang,glb,arq.fdis);
        %arq.ftes=esff'./area';
        %arq.fen=F'*arq.fdis;
        arq.flmbd=Lambd(1);
        arq.ilmbd=Lambd(1);
    end
    
    %VIBRAÇÃO
    if arq.modos(3)==0
        
        arq.ffreq=0;
        arq.ifreq=0;
        arq.escvib=0;
    else
        disp ('Freq')
        freq=freq(arq.modos(3));
        arq.escvib=1;
        arq.D=du(:,arq.modos(3));
        disp(freq)
        %arq.x=area;%/Lambd;
        %load trussdata.mat
        %esff = tresf(arq.x,comp,els,ang,glb,arq.fdis);
        %arq.ftes=esff'./area';
        %arq.fen=F'*arq.fdis;
        arq.ffreq=freq(1);
        arq.ifreq=freq(1);
    end

    
    arq.Resumo.funcCount=1;
    arq.ites=arq.ites';
    arq.ftes=arq.ites;
    arq.fvol=arq.ivol;
    arq.fen=arq.ien;

%     disp('Deslocamentos')
%     disp(arq.fdis(:)')
%     disp('Tensoes')
%     disp(arq.ftes(:)')
%     disp('Areas')
%     disp(arq.x0);
%    disp(arq)
end

%  
%  [arq.fdis,arq.ftes,esf,arq.fvol,alph,fn,kn,kni,Z]=fesolrm2(arq.x,SN,N);
disp('tempo')
disp(arq.t)

end