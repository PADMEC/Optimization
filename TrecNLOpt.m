%function Trec2OptrussOt
disp(' Otimização Estrutural ')

% Programa para analise de treliças planas pelo metodo dos deslocamentos

        global  reliaboptimiz failure_sets restart_parm
        global  HYPLAS_flag command_string
        global  modos POD_Z  xol Xs P0 para itp
        Xs = [];
restart_parm=[0,0,0];
if isempty(HYPLAS_flag);
    HYPLAS_flag = 0;
end
% Entrada de Dados relacionado com a Geometria dos elementos
%
% forneça as coordenadas dos nos
% load  C:\MATLAB7\work\arq0

T11=cputime;
if HYPLAS_flag
    [hyplas_command] = Hyplas_analysis( arq.x0, hyplasprojname,0,RandVarDist(2,:));
    [Fs,Us,Sigs,ien,vol  ] = Hyplas_analysis( arq.x0, hyplas_command, -3, RandVarDist(2,:));
%     [feFs,feUs,feSigs,feien,feivol ] = Hyplas_analysis( arq.x0, hyplas_command,-3, RandVarDist(2,:));
    varginame={arq.x0, hyplasprojname,0,0};
    if RBtype==0
        RBtype=3;
    end
    command_string=hyplas_command;
else
    [Fs,Us,Sigs,ien,vol,T,floc,u_flamb,Lambd,D,freq]=NLfesol(area,props,fext,glb,tipo);
end
T2=cputime;
Tfem=T2-T11
arq.vol0=vol;

%figure, plot_sig_vet(arq.x0(arq.arex),fext,Us(:,end),Sigs(:,end),1)
% Lambd
%figure, plot_sig_vet(arq.x0(arq.arex),fext,u_flamb,Sigs(:,end),-1)
global convergence_RBM i_convR f_convR
convergence_RBM=0;
if RBtype==1||RBtype==2
    %[D,Fx,Toff,x0s,k_iter]=OfflinePOD(arq,props,fext,glb,tipo);
    N=noffT; %n of non-lin offline analyses 
    if HYPLAS_flag
        [Dt,Fx,Cs,Toff,x0s,k_iter]=HplsOfflinePOD(N,arq,funame,varginame);
    else
        [Dt,Fx,Zkn,Cs,Toff,x0s,k_iter]=OfflinePOD(arq,props,fext,glb,tipo,N);
    end
    
    if ioffT<noffT
        Convergence_ReducModel
    end
    
    if RBtype==1 %POD
        lamb_tol=1e-5;
%         Ps=D;
%       P0=mean(D,2); Ps = D - P0*ones(1,size(D,2));
        [POD_Z,nPZ,V]=PODfunction(Dt,lamb_tol);
    else % RBM
        POD_Z = D;
    end
end % end RBM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pre-proc
itera = 0;
perturb=1e-6;
tpres=arq.tpres;
ndvab=length(arq.x0);
vlb =  arq.xmin* ones(1,ndvab);
vub =  arq.xmax* ones(1,ndvab);

if HYPLAS_flag
    write_Zfile(POD_Z,'RB_file.dat')
    vargs={x, hyplas_command,-RBtype,RandVarDist(2,:)};
%   [ pFs, pUs, pSigs, pien,pivol  ] = Hyplas_analysis( vargs{:});
else
    vargs={area,props,fext,glb,tipo,RBtype,POD_Z};
end

global AVeq BVeq
AVeq=[];BVeq=[];
if ~HYPLAS_flag

% arq.ivol=sum(area.*comp.*ro);
    %[ndvab,link,lpdva,perturb]=sendat;
     link=[arq.arex(:) ones(nel,1)];
     
    %%%%%%%%%%%%%%%%%
    % Volume constraint
    if arq.tpres(1)
        AVeq=zeros(1,ndvab);
        for iv=1:ndvab
            AVeq(1,iv)=sum(props{2}(link(:,1)==iv));
        end
        BVeq=AVeq*arq.x0(:);
        arq.tpres(1)=0;
    end
    %%%%%%%%%%%%%%%%%
    
%   entrada de Dados relacionado com a otimização
 		%[tobj,tpres,vlb,vub,x0,clb,cub]=optdat(area,lpdva,ndvab,nno);
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
    varginame={arq.x0, hyplas_command,-RBtype,RandVarDist(2,:)};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stochastic Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Stochast_on
    ft = RandVarDist(1,:);
    Mu = RandVarDist(2,:);
    Su = RandVarDist(3,:);
    Cov= Correl.*(Su'*Su);
    N1 = DinamcVar;
    
    if Stochast_moment||strcmp(Method,'R2DO')
        %% Statistics Computation
        ntest=size(X_ops,1);
        for iopt=1:ntest
            %xopt=arq.x0(iopt,:);
            xopt=X_ops(iopt,:)
            if HYPLAS_flag
                vargs{1}=xopt;
            else
                area = xopt(arq.arex(:,1));
                vargs{1}=area;
            end
            tic
            [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc( Mu,Cov,ft,S_fun,Method,N1,vargs{:});
            toc
            %v_out = {vol,u,vonmis,sn,dvol,du,dsig,den};
            if ntest>1
                fun='myNLfun';
                funobj(iopt,:)=feval(fun,xopt,vargs{:});
                fun='myNLcon';
                funconst=feval(fun,xopt,vargs{:});
                [const_iopt(iopt),iactive(iopt)]=max(funconst);
            end
        end
        % Test of Random Var Importance
        if test_importance
            Importance_factors
            pause
        end
    end
    
    if ~Stochast_moment||strcmp(Method,'R2DO')
        %% Reliability Analysis
        reliaboptimiz = 'RIA';
%          reliaboptimiz = 'PMA';
        
        if ~Stochast_moment
            N1 = 2;
            [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = ...
                StochastiCalc( Mu,Cov,ft,S_fun,'PC',N1,vargs{:});
        end
        %Approximated
        [beta,dbeta,failure_sets] = Reliab_approx(reliaboptimiz,...
            M_out{G_const},dM_out{G_const},S_out{G_const},dS_out{G_const});
%         failure_sets=14;
        if ~isempty(failure_sets)
            if strcmp(reliaboptimiz,'PMA')
                [Vopt,betasi,dbetasi,Converg]=PMA_RBDO(Mu,Cov,ft,G_fun,vargs{:});
            else
                [PF,beta, dbeta, v_outPP] = Reliab...
                    (Mu,Cov,ft,G_fun,Method,DinamcVar,1,vargs{:});
            end
        end
        %disp([xvalu, dbeta'])
    end
end




if arq.tobj==0      %  MULTI-OBIJETIVO
    %% Vector Optimization
    
    MOptTypes={'ws','mmx','nbi','nnc','nbim','nncm'};
    imethods=1:4;%4;%2;%1:4;%[3,4,2,1];
    ptsp=10;
    for i=imethods
        MOptype=MOptTypes{i}
        %xol = xpdva;
        t1 = tic;
        %[x, f, Cnvrg, fcount,Tot_Time]=Mfmincon([1 1 1], [1 1 1]*(-9999),[1 1 1]*(9999), 'fun33mo', 'con33mo', 10, MOptype);
        [x, f, Cnvrg, fcount,Tot_Time]=Mfmincon(xpdva, vlb,vub, 'myNLfun', 'myNLcon', ptsp, MOptype,xpdva,fext,glb,link);
        %arq=moptruss(arq,plt, arq.x0,vlb,vub,xpdva,'MObj',1,fext,glb,link);
        %                                    ,{  varargin  }
        times_all(i)=toc(t1);
        %load  C:\MATLAB7\work\arq0
        %Save method results
        Resum{i}={f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp,mo};
        arq.t(i) = Tot_Time;
        MOptype,f, Cnvrg,times_all(i)
    end
    save(['ResumMO_' Method '_' num2str(ceil(cputime))],'Resum','imethods','mo','tpres','times_all')
    arq.x = x(ceil(ptsp/2),arq.arex);
    times_all
    nobj=size(Resum{imethods(1)}{1},2);
    figure, hold on
    for i=imethods
        switch i
            case 1, pt='s' ;leg{i}='WS';
            case 2, pt='^g';leg{i}='Min Max';
            case 3, pt='ok';leg{i}='NBI';
            case 4, pt='pr';leg{i}='NNC';
            case 5, pt='.k';leg{i}='NBIm';
            case 6, pt='.r';leg{i}='NNCm';
        end
        if nobj==2
            plot(Resum{i}{1}(:,1),Resum{i}{1}(:,2),pt)
        elseif nobj==3
            grid on
            plot3(Resum{i}{1}(:,1),Resum{i}{1}(:,2),Resum{i}{1}(:,3),pt)
        end
        %[TTs(i),TFCs(i),eveness(i),Par_area1(i),Par_area(i),nnonP(i),ndpP(i)]=post_proc(Resum{i}{:});
        
    end
    legend(leg{imethods})
%     % ROBUST OPTIMIZATION
%     xlabel('Mean (N/cm²)')
%     ylabel('Standard Deviation (N/cm²)')

    %for i=imethods
        [TTs,TFCs,eveness,Par_area1,Par_area,nnonP,ndpP]=post_proc(Resum);
    %end
    %PNG_n_CLOSE
    Resume_Tablet=...
    [TTs/60;TFCs;nnonP;eveness;Par_area];
    f=Resum{3}{1};
    % Plot Pareto Designs
    for i=1:ptsp
        xi=x(i,arq.arex);
        disp([i,x(i,:)])

        if HYPLAS_flag
            vargs={xi, hyplas_command,-RBtype,RandVarDist(2,:)};
            [Fs,Us,Sigs,arq.en,arq.vol  ] = Hyplas_analysis(vargs{:});

        else
            vargs={xi,props,fext,glb,-tipo,RBtype,POD_Z};
            % FEM (fully)
            [Fs,Us,Sigs,arq.en,arq.vol,T]=NLfesol(xi,props,fext,glb,tipo);
            % RBM - POD
            [Fs,Us,Sigs,arq.en,arq.vol,T]=NLfesol(vargs{:});
        end
        Maxu(i) = max(abs(Us(:,end)));
        Maxs(i) = max(abs(Sigs(:,end)));
        
        if Stochast_on
 
            if Stochast_moment||strcmp(Method,'R2DO')
                [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc...
                    ( Mu,Cov,ft,S_fun,Method,N1,vargs{:});
            end
            M=M_out{G_const}(set_Failure_FORM);
            S=S_out{G_const}(set_Failure_FORM);
            btsa(i)=max(G_limits-abs(M))/S;

            if ~Stochast_moment||strcmp(Method,'R2DO')
                if strcmp(reliaboptimiz,'PMA')
                    [Vopt,beta,dbeta,Converg]=PMA_RBDO(Mu,Cov,ft,G_fun,vargs{:});
                else
                    [PF,beta, dbeta, v_outPP] = Reliab...
                        (Mu,Cov,ft,G_fun,Method,DinamcVar,1,vargs{:});
                end
                bts(i)=beta;
                PFs(i)=PF;
            end
        end
    end
    figure, plot(Maxu([2:10,1]),Maxs([2:10,1]),'.-')
    xlabel('Max u')
    ylabel('Max sig')
    if ~Stochast_moment||strcmp(Method,'R2DO')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PLOT Pareto x Beta index
        figure,hold on
        [AX,H1,H2]=plotyy(f([2:10,1],1),f([2:10,1],2),f([2:10,1],1),bts([2:10,1]),'plot');
        set(AX,'FontWeight','bold','FontSize',14);
        set([H1,H2],'Marker','.','MarkerSize',15,'LineWidth',2);
        
        % beta aproximated (G-M)/S
        hold(AX(2),'on')
        plot(AX(2),f([2:10,1],1),btsa([2:10,1]),'o-','MarkerSize',10,'LineWidth',2)
        
        % Labels and Legends
        HL=legend('Pareto','{\it{\beta}} FORM','{\it{\beta}} Aprox');
        HLb(1)=xlabel(AX(1),'{\it{F}}_{1}');
        HLb(2)=ylabel(AX(1),'{\it{F}}_{2}');
        HLb(3)=ylabel(AX(2),'{\it{\beta}}');
        set([HL;HLb(:)],'FontWeight','bold','FontSize',14)
        grid on
%         hold(AX(1),'on')
%         plot(AX(1),ff([2:10,1],1),ff([2:10,1],2),'x-','MarkerSize',10,'LineWidth',2)
%         Title('\beta','FontWeight','bold','FontSize',14)
    end
    
    
else
    %% Scalar Optimization

    %   Uso da biblioteca de otimização (fmincon)
    %
   options=optimset;
   %options=optimset(options,'LargeScale','on');'SubproblemAlgorithm','cg'
   options=optimset(options,'GradObj','on','GradConstr','on','MaxFunEvals',40000,'Algorithm' ,'sqp');%'active-set');%'interior-point');%,
   options=optimset(options,'MaxIter',10000,'Display','iter','DerivativeCheck','off', 'LargeScale','on');%,'MaxTime',9);
   %options=optimset(options,'LargeScale','off','Algorithm','interior-point');
 
    %       
    %	     Solução final
    %
    disp('Otimizando...')
    tic
    [xopt,ffobj,arq.Converg,arq.Resumo]=fmincon('myNLfun',arq.x0,[],[],AVeq,BVeq,vlb,vub,'myNLcon',options,xpdva,fext,glb,link);
    arq.t = toc;

    xopt
    ffobj
    arq.x=xopt(arq.arex);
end
%[arq.fdis,arq.ftes,arq.fen,arq.fvol,floc,du,Lambd,D,freq]=fesol(arq.x,props,fext,glb,[1 1 0]);

if HYPLAS_flag
    % VARIATION PARAMETRIC PLOT
%     x1=100:10:477;
%     for i=1:length(x1)
%         arq.x=[x1(i),100]
%     vargs={arq.x, hyplas_command,-RBtype,RandVarDist(2,:)};
%     [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc...
%                         ( Mu,Cov,ft,S_fun,Method,N1,vargs{:});
%         Sis(i)=M_out{3}(264,end)
%         Uis(i)=M_out{2}(3,end);
%         SSis(i)=S_out{3}(264,end)
%         SUis(i)=S_out{2}(3,end);
% 
%         [mS(i),kms(i)]=max(M_out{3}(:,end));
%         [mU(i),kmu(i)]=max(M_out{2}(:,end));
%         [smS(i),ksms(i)]=max(S_out{3}(:,end));
%         [smU(i),ksmu(i)]=max(S_out{2}(:,end));
%         [mG(i),kmg(i)]=max(M_out{2}(:,end)+3.2*S_out{2}(:,end));
%         figure,plot(Sis,SSis,'.')

    [Fs,Us,Sigs,arq.en,arq.vol  ] = Hyplas_analysis(arq.x, hyplas_command,-RBtype,RandVarDist(2,:));
%         Sis(i)=Sigs(264,end);
%         Uis(i)=Sigs(3,end);
%         [mS(i),kms(i)]=max(Sigs(:,end));
%         [mU(i),kmu(i)]=max(Us(:,end));
%     end
else
    [Fs,Us,Sigs,arq.ien,arq.ivol,T,floc,u_flamb,Lambd,D,freq]=NLfesol(arq.x,props,fext,glb,tipo);
end

if Stochast_moment||strcmp(Method,'R2DO')
    vargs={arq.x,props,fext,glb,tipo,RBtype,POD_Z};
    [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc...
        ( Mu,Cov,ft,S_fun,Method,N1,vargs{:});
            en = v_out0{4};
            Men = M_out{4};
            Sen = S_out{4};
            [Men,Sen]
end

%% PLOT STRUCTURE
if HYPLAS_flag
%     varginame={arq.x, hyplasprojname,0,RandVarDist(2,:)};
    plot_HYPLAS(funame,{arq.x, hyplasprojname,0,RandVarDist(2,:)})
else
    figure, plot_sig_vet(arq.x,fext,Us(:,end),Sigs(:,end),1)
end

arq.Resumo
arq.Converg

disp('tempo de otimizacao')
disp(arq.t)
arq.t=sum(arq.t);

% save  C:\MATLAB7\work\arq0 arq plt

