%function Trec2OptrussOt
disp('Strutural Optimization via HYPLAS')
 
% Programa para analise de treliças planas pelo metodo dos deslocamentos

        global modos POD_Z
        global  para con
        global   xol Xs P0
        global itp command_string
        Xs = [];
%
% Entrada de Dados relacionado com a Geometria dos elementos
%
% forneça as coordenadas dos nos
% load  C:\MATLAB7\work\arq0
if ~Stochast_on
    RandVarDist=zeros(3,0);
end
[hyplas_command] = Hyplas_analysis( arq.x0, hyplasprojname,0,RandVarDist(2,:));

T11=cputime;
%[Fs,Us,Sigs,arq.ien,arq.ivol,T]=NLfesol(area,props,fext,glb,tipo);
[Fs,Us,Sigs,ien,vol  ] = Hyplas_analysis( arq.x0, hyplas_command, -3, RandVarDist(2,:));
T2=cputime;
Tfem=T2-T11
arq.vol0=vol;
command_string=hyplas_command;

funame = 'Hyplas_analysis';
varginame={arq.x0, hyplasprojname,0,0};

if RBtype==0
    RBtype=3;
end

if RBtype==1||RBtype==2
    %[D,Fx,Toff,x0s,k_iter]=OfflinePOD(arq,props,fext,glb,tipo); 
    noffT = 10;
    ioffT = 10;
    for in = ioffT:noffT % BASIS CONVERGENCE
        nOff=in; %n of non-lin analyses offline
        [D,Fx,Cs,Toff,x0s,k_iter]=HplsOfflinePOD(nOff,arq,funame,varginame);
        
        if RBtype==1
            lamb_tol=1e-6;
            P0=mean(D,2); Ps = D - P0*ones(1,size(D,2));
            [POD_Z,nPZ,V]=PODfunction(Ps,lamb_tol);
        else
            POD_Z = D;
        end
        
        write_Zfile(POD_Z,'RB_file.dat')
        T11=cputime;
        [Fsn,Usn,Sigsn,ienn,voln  ] = Hyplas_analysis( arq.x0, hyplas_command,-RBtype, RandVarDist(2,:));
        T2=cputime;
        Trbm=T2-T11
        %figure, hold on, plot(Sigs(:,end),'o'), plot(Sigsn(:,end),'.r')

        %area=area/4;
        %[pFs,pUs,pSigs,pien,pivol,pT]=NLfesol(area,props,fext,glb,tipo,RBtype,POD_Z);
        
        %Test POD
%         [feFs,feUs,feSigs,feien,feivol,feT]=NLfesol(area,props,fext,glb,tipo,0);
%         [rFs,rUs,rSigs,rien,rivol,rT]=NLfesol(area,props,fext,glb,tipo,RBtype,D);
%         [mr,ir]=max(abs((rUs-feUs)));disp(norm(rUs-feUs)/norm(feUs))
%         [mp,ip]=max(abs((pUs-feUs)));disp(norm(pUs-feUs)/norm(feUs))
%         
%         figure, plot(rUs-feUs,'.-'), hold on, plot(pUs-feUs,'o-k'), plot(feUs,'xr')

        if ioffT<noffT
            errP(in-ioffT+1)=arq.ien-pien
            errR(in-ioffT+1)=arq.ien-rien
            errPd(in-ioffT+1)=max(abs(pUs-Us));
            errRd(in-ioffT+1)=max(abs(rUs-Us));
            nPz(in-ioffT+1)=size(POD_Z,2);
        end
    end
    if ioffT<noffT
        figure, hold on
        plot([ioffT:noffT],errP)
        plot([ioffT:noffT],errR)
        title ('Con Error')
        xlabel('N')

        figure, hold on
        plot([ioffT:noffT],errPd)
        plot([ioffT:noffT],errRd)
        title ('Displ Error')
        xlabel('N')
    end
%vargs={area,props,fext,glb,tipo,RBtype,POD_Z};
end

% varginame={x, hyplasprojname,RBtype,RandVarDist(2,:)};
% plot_HYPLAS(funame,varginame)
varginame={x, hyplas_command,-RBtype,RandVarDist(2,:)};

%Testing optimum points
% varginame={xopt(2,:), hyplas_command,-RBtype,RandVarDist(2,:)};
if Stochast_on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stochastic

    ft = RandVarDist(1,:);
    Mu = RandVarDist(2,:);
    Su = RandVarDist(3,:);
    Cov= Correl.*(Su'*Su);
%     G_fun = 'HyplasLimitState';
%     S_fun = 'Stoch_HyplasAnalysis';
    N1 = DinamcVar;
%         if strcmp(Method,'FORM')
%                 DinamcVar=i_contin;
%         end
    if Stochast_moment||strcmp(Method,'R2DO')||strcmp(Method,'R2DOmc')
        % Statistics Computation
        xopt=x;
        ntest=size(xopt,1);
        for iopt=1:ntest
%             xopt=arq.x0(iopt,:)
            varginame{1}=xopt(iopt,:);

            T11=cputime;
            [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc( Mu,Cov,ft,S_fun,Method,N1,varginame{:});
            T2=cputime;
            Tstochmom=T2-T11
        %v_out = {vol,u,vonmis,sn,dvol,du,dsig,den};
        
%         dS_out{4}',Sen=S_out{4}
%         area0=area;
%         for i=1:length(area)
%             area=area0;
%             area(i)=area(i)+0.5;
%             vargs={area,props,fext,glb,tipo,RBtype,POD_Z};
%             [ iv_out0, iM_out, iS_out, idv_out0, idM_out, idS_out ] = StochastiCalc( Mu,Cov,ft,S_fun,Method,N1,varginame{:});
%             newSen(i)=iS_out{4};
%             newMen(i)=iM_out{4};
%         end
        
            Men = M_out{4};
            Sen = S_out{4};
            Riopt(iopt,:)=[Men,Sen];
            Riopt2(iopt,:)=[max(abs(M_out{2})), max(abs(S_out{2}))];
            const_iopt(iopt)=max(abs(M_out{3})+3*S_out{3});
        end
        % Test of Random Var Importance
        if test_importance
            stds=Su;
%             designRV=RandVarProp(1,:)==1;
%             componts={rv_data_comp{designRV}};
%             stds(designRV)=Su(designRV).*arq.x0([componts{:}]);

            en = v_out0{4};
            Men = M_out{4};
            Sen = S_out{4};

            den = dv_out0{4};
            dMen = dM_out{4};
            dSen = dS_out{4};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % STRESS IMPORTANCE FACTOR COMPUTATION
            i_IFC=2;
            [VM,imax] = max(abs(M_out{i_IFC}));
%             [VM,imax] =max(abs(M_out{i_IFC})+3*S_out{i_IFC})
            en = v_out0{i_IFC}(imax);
            Men = M_out{i_IFC}(imax);
            Sen = S_out{i_IFC}(imax);
            
            den = dv_out0{i_IFC}(:,imax);
            dMen = dM_out{i_IFC}(:,imax);
            dSen = dS_out{i_IFC}(:,imax);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            DeltEnergy=[dMen'; dSen'].*([1/Men;1/Sen]*stds);
            DeltEnergy=[den';dMen'; dSen'].*([1/en;1/Men;1/Sen]*stds);
            figure,bar(abs(DeltEnergy'))
            legend('Deterministic','Mean','S.D.')
            title(['R.V. importance for \it{x} = ' num2str(arq.x0,'%0.5g ')])
            xlabel('Random Variables')
            ylabel('Output Variation (%)')
            disp('paused...'),pause
        end
        
    end
    
    if ~Stochast_moment||strcmp(Method,'R2DO')||strcmp(Method,'R2DOmc')
        % Reliability Analysis
%         [PF,beta, dbeta,Fs,Us,Sigs,en,vol,dF,du,dsig,den,dvol] =
%         Reliab(Mu,Cov,ft,G_fun,Method,DinamcVar,1,varginame{:});
        T11=cputime;
        [PF,beta,dbeta,varargout] = Reliab(Mu,Cov,ft,G_fun,Method,DinamcVar,1,varginame{:});
        T2=cputime;
        Treliab=T2-T11
        %disp([xvalu, dbeta'])
    end
end



% 
% [arq.idis,arq.ites,arq.ien,vol,floc,du,Lambd,D,freq]=fesol(area,props,fext,glb,[1 1 0]);
% arq.ilmbd=Lambd(1);
% arq.ivol=sum(area.*comp.*ro);
%
% Execute a opção desejada:
%		3--- OTIMIZAÇ+O
itera = 0;
perturb=1e-6;
tpres=arq.tpres;
ndvab=length(arq.x0);
vlb =  arq.xmin* ones(1,ndvab);
vub =  arq.xmax* ones(1,ndvab);

if ~HYPLAS_flag
%
		%[ndvab,link,lpdva,perturb]=sendat;
         link=[arq.arex(:) ones(nel,1)];
%   entrada de Dados relacionado com a otimização
%
 		%[tobj,tpres,vlb,vub,x0,clb,cub]=optdat(area,lpdva,ndvab,nno);
%
        %arq.tobj=2;%    1---VOLUME  2----ENERGIA   3--Desl.Tot
        %4--Desl.Maxim   5--Desl.Esp   6--Tes.Max
        %arq.tpres=1;%    1---PESO    2----TENSAO     3---Tens e Desloc
        %4---Desloc.
        
% [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc...
%                     ( Mu,Cov,ft,S_fun,Method,N1,vargs{:});
                
    disp('energy mean,std:')
    M_out{4},S_out{4}
        
        fs=4; %     Func. Subs. - 1. Soma Pond - 2. Min-Max - 3.NBI - 4. NNC
        %pareto=1;   % Varia o beta e acha 
         
 
%
%   Especfc
%

    [fde,fte,tpres]=fobrest(arq,length(fext),area,arq.E,arq.ro,props{2},[],[],[]);
    para = [perturb;ndvab;arq.tobj;tpres(1);itera;1;1];
    % props = [ang;comp;els;ro];
    modos=arq.modos;

    xpdva =  arq.x0;
    arq.t=[];
    xol = xpdva*0;
else
    itp=find(tpres);
    xpdva = arq.x0;
    fext=[];
    glb=[];
    link=1:ndvab;
    para = [perturb;ndvab;arq.tobj;tpres(1);itera;1;1];
    nn=size(Sigs,1);
    [fde,fte,tpres]=fobrest(arq,size(Us,1),0,0,0,0);
end

varginame={arq.x0, hyplas_command,-RBtype,RandVarDist(2,:)};
if arq.tobj==0      %  MULTI-OBIJETIVO
    MOptTypes={'ws','mmx','nbi','nnc','nbim','nncm'};
    imethods=1:4
    ptsp=10;
    for i=imethods
        MOptype=MOptTypes{i}
        %xol = xpdva;
        [x, f, Cnvrg, fcount,Tot_Time]=Mfmincon(xpdva, vlb,vub, 'myNLfun', 'myNLcon', ptsp, MOptype,xpdva,fext,glb,link);
        %arq=moptruss(arq,plt, arq.x0,vlb,vub,xpdva,'MObj',1,fext,glb,link);
        %                                    ,{  varargin  }
        %load  C:\MATLAB7\work\arq0
        %Save method results
        Resum{i}={f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp};
        arq.t(i) = Tot_Time;
        MOptype,f, Cnvrg
    end
    save(['ResumMO_' Method '_' num2str(ceil(cputime))],'Resum','imethods','mo','tpres')
    arq.x = x(ceil(ptsp/2),arq.arex);
    figure, hold on
    for i=imethods
        switch i
            case 1, pt='*' ;leg{i}='WS';
            case 2, pt='^g';leg{i}='Min Max';
            case 3, pt='ok';leg{i}='NBI';
            case 4, pt='+r';leg{i}='NNC';
        end
        plot(Resum{i}{1}(:,1),Resum{i}{1}(:,2),pt)
        %[TTs(i),TFCs(i),eveness(i),Par_area1(i),Par_area(i),nnonP(i),ndpP(i)]=post_proc(Resum{i}{:});
        
    end
    legend(leg{imethods})
    % ROBUST OPTIMIZATION
    xlabel('Mean')
    ylabel('Standard Deviation')

    for i=imethods
        [TTs(i),TFCs(i),eveness(i),Par_area1(i),Par_area(i),nnonP(i),ndpP(i)]=post_proc(Resum{i}{:});
    end
    %PNG_n_CLOSE
    Resume_Tablet=...
    [TTs/60;TFCs;nnonP;eveness;Par_area];

     varginame={x(1,:), hyplasprojname,RBtype,RandVarDist(2,:)};
     plot_HYPLAS(funame,varginame)
     varginame={x(2,:), hyplasprojname,RBtype,RandVarDist(2,:)};
     plot_HYPLAS(funame,varginame)


%      plot(Resum{2}{1}(:,1),Resum{2}{1}(:,2),'^g')
%     hold on, plot(Resum{3}{1}(:,1),Resum{3}{1}(:,2),'ok')
%     hold on, plot(Resum{4}{1}(:,1),Resum{4}{1}(:,2),'+r')

    for i=1:ptsp
        xi=x(i,arq.arex);
        varginame{1}=xi;
        % FEM (fully)
%         [Fs,Us,Sigs,arq.en,arq.vol,T]=NLfesol(xi,props,fext,glb,tipo);
        % RBM - POD
%         [Fs,Us,Sigs,arq.en,arq.vol,T]=NLfesol(vargs{:});
%         [Fs,Us,Sigs,arq.en,arq.vol  ] = Hyplas_analysis( xi, hyplas_command,-3,RandVarDist(2,:));
        [Fs,Us,Sigs,arq.en,arq.vol  ] = Hyplas_analysis(varginame{:});



        Maxu(i) = max(abs(Us(:,end)));
        Maxs(i) = max(abs(Sigs(:,end)));
        if Stochast_on
 
%             if Stochast_moment||strcmp(Method,'R2DO')
%                 [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc...
%                     ( Mu,Cov,ft,S_fun,Method,N1,varginame{:});
% 
%             end

            if ~Stochast_moment||strcmp(Method,'R2DO')
                [PF,beta, dbeta,Fs,Us,Sigs,en,vol,T,floc,u_flamb,lamb,D,freq,dF,du,dsig,den,dvol] = Reliab...
                    (Mu,Cov,ft,G_fun,Method,DinamcVar,1,varginame{:});

                bts(i)=beta;
                PFs(i)=PF;
            end
        end
    end
    figure, plot(Maxu,Maxs,'.')
    xlabel('u')
    ylabel('Seff')
else
    % Scalar Optimization (Uni-objective problem))
    %   Uso da biblioteca de otimização (fmincon)
    %            		  
   options=optimset;
   %options=optimset(options,'LargeScale','on');'SubproblemAlgorithm','cg'
   options=optimset(options,'GradObj','on','GradConstr','on','MaxFunEvals',40000,'Algorithm' ,'sqp');%'active-set');%'interior-point');%'sqp');%,
   options=optimset(options,'MaxIter',10000,'Display','iter','DerivativeCheck','off', 'LargeScale','off');%,'MaxTime',9);
   %options=optimset(options,'LargeScale','off','Algorithm','interior-point');
 
    %       
    %	     Solução final
    %
    disp('Otimizando...')
    tic
    [xopt,ffobj,arq.Converg,arq.Resumo]=fmincon('myNLfun',arq.x0,[],[],[],[],vlb,vub,'myNLcon',options,xpdva,fext,glb,link);
    arq.t = toc;

    xopt
    ffobj
    x=xopt;
%     MainPostproc
%     checkbox = zeros(1,11);
%     checkbox([2,4])=1;
%     PlotResults(checkbox)
    arq.x=xopt(arq.arex);
end

[Fs,Us,Sigs,ien,vol  ] = Hyplas_analysis( arq.x, hyplasprojname,0,RandVarDist(2,:));
[RBFs,RBUs,RBSigs,RBien,RBvol  ] = Hyplas_analysis( arq.x, hyplasprojname,RBtype,RandVarDist(2,:));

varginame={arq.x, hyplasprojname,0,RandVarDist(2,:)};
plot_HYPLAS(funame,varginame)

%[arq.fdis,arq.ftes,arq.fen,arq.fvol,floc,du,Lambd,D,freq]=fesol(arq.x,props,fext,glb,[1 1 0]);
%[Fs,Us,Sigs,arq.en,arq.vol,T]=NLfesol(arq.x,props,fext,glb,RBtype);
%{   u,      sig,  en,  vol,   floc,du,Lambd,uf,sigf,esff}
%arq.fdis=uf;
%     [dis,tes,en,vol,floc,du,Lambd,uf,sigf,esff]=fesol(arq.x,props,fext,glb,1);
%     arq.fdis=uf;
%     arq.ftes=sigf;
%disp('Lambda = ');disp(Lambd(1))%:ceil(nglb/2))'
%arq.flmbd=Lambd(1);
if Stochast_moment||strcmp(Method,'R2DO')
    varginame={arq.x, hyplas_command,-RBtype,RandVarDist(2,:)};
%     vargs={arq.x,props,fext,glb,tipo,RBtype,POD_Z};
    [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc...
        ( Mu,Cov,ft,S_fun,Method,N1,varginame{:});
            en = v_out0{4};
            Men = M_out{4};
            Sen = S_out{4};
            [Men,Sen]
end

    %figure, plot_sig_vet(arq.x,fext,Us(:,end),Sigs(:,end),1)

arq.Resumo
arq.Converg

disp('tempo de otimizacao')
disp(arq.t)
arq.t=sum(arq.t);

% save  C:\MATLAB7\work\arq0 arq plt

