function [fob,gob,fre,gre]=optNLsol(props,fext,glb,link,tobj,tres,xvalu,ndvab,~,~)
    % Analysis
    global modos RBtype POD_Z Method
    
    % Stochastic Variables
    global Stochast_on RandVarDist Correl Stochast_moment DinamcVar 
    
    % Reliabirelity Analysys Variables
    global G_const failure_sets betatarg reliaboptimiz set_Failure_FORM

    %HYPLAS Variables
    global HYPLAS_flag command_string 
    
%   Multi-objective case (test nbi)
%%%   NBI   %%%
global t fs
if tobj == 0
    if isempty(fs)
        fs=0;
    end
    %global mo
    %fmof=mo(1,3);
    if fs==3||fs==5
        %%%   NBI   %%%
        t=xvalu(end);
        xvalu(end)=[];
    end
end
    
if ~HYPLAS_flag
    
    % Fornece valores para as funcões objetivo e restrição e suas derivadas
    %
%     ang  = props(1,:);
%     comp = props(2,:);
%     els  = props(3,:);
%     ro  = props(4,:);
%     X = xvalu;

    %area pre-calculation
    area = xvalu(link(:,1)).*link(:,2)';

    %input of the nonlin analysis
    vargs={area,props,fext,glb,-modos,RBtype,POD_Z};
    G_fun = 'StructLimitState';
    S_fun = 'Stoch_Analysis';
else
    vargs={xvalu, command_string,-RBtype};
    G_fun = 'HyplasLimitState';
    S_fun = 'Stoch_HyplasAnalysis';
end

if Stochast_on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stochastic Analysis
    ft = RandVarDist(1,:);
    Mu = RandVarDist(2,:);
    Su = RandVarDist(3,:);
    Cov= Correl.*(Su'*Su);
    N1 = DinamcVar;
%         if strcmp(Method,'FORM')
%                 DinamcVar=i_contin;
%         %else
%         end

    if Stochast_moment||strcmp(Method,'R2DO')||strcmp(Method,'R2DOmc')
        %% Compute Stochastic Moments
        [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc...
            ( Mu,Cov,ft,S_fun,Method,N1,vargs{:});
        %{vol,u,vonmis,sn,dvol,du',dsig',den'};
        %v_out = {vol,u,vonmis,sn,dvol,du,dsig,den};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% test plot
        %sig=M_out{3};u=M_out{2};
        %tp=1; title_str='SigM';
        %figure, plot_sig_vet(area,fext,u*50,sig,tp)
        %title(title_str)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         dx=1e-5;iv=1;kk=1;iout=3;
%         xvalu(iv)=xvalu(iv)+dx;
%         area = xvalu(link(:,1)).*link(:,2)';
%         vargs={area,props,fext,glb,-modos,RBtype,POD_Z};
%         [ v_out01, M_out1, S_out1, dv_out01, dM_out1, dS_out1 ] = StochastiCalc...
%             ( Mu,Cov,ft,S_fun,Method,N1,vargs{:});
%         [(M_out1{iout}(kk)-M_out{iout}(kk))/dx,dM_out{iout}(iv,kk)]

        
    end
    if ~Stochast_moment||strcmp(Method,'R2DO')||strcmp(Method,'R2DOmc')
        %% Reliability Based Design Optimization
        
        %Approximated
        [beta,dbeta,failure_sets] = Reliab_approx(reliaboptimiz,M_out{G_const},dM_out{G_const},S_out{G_const},dS_out{G_const});
        
        if ~isempty(failure_sets)
            if strcmp(reliaboptimiz,'PMA')
                % PMA
                [Vopt,betasi,dbetasi,Converg]=PMA_RBDO...
                    (Mu,Cov,ft,G_fun,vargs{:});
                
                if ~isempty(set_Failure_FORM)
                    % Set a specific failure component
                     beta =  betasi';
                    dbeta = dbetasi;
                else
                    % Set all failures func
                     beta(failure_sets) = betasi';
                    dbeta(:,failure_sets) =dbetasi;
                end
                 ReliabC =  beta;
                dReliabC = dbeta;
            else
                % RIA
                [~,betasi,dbetasi] = Reliab...
                    (Mu,Cov,ft,G_fun,Method,DinamcVar,1,vargs{:});

                if ~isempty(set_Failure_FORM)
                    % Set a specific failure component
                     beta =  betasi';
                    dbeta = dbetasi;
                else
                     beta(failure_sets) = betasi;
                    dbeta(:,failure_sets) =dbetasi;
                end
                 ReliabC = (betatarg - beta);
                dReliabC = -dbeta;

            end
        else
%             betas(betas==0)=inf;
            % Compute failure probabilities
            Ps = normcdf(-betas );
            Ps(isnan(Ps))=0;
            % Compute union of failure probabilities
            PF = sum(Ps);
            if (PF<1e-16)
                [beta,kbm] = min(betas);
                dbeta=dbetas(:,kbm);
            elseif(PF>(1-1e-16))
                [beta,kbm] = max(betas);
                dbeta=dbetas(:,kbm);
            else
                beta = - inv_norm_cdf(PF);

                % Sens
                ndv=length(xvalu);
                dPis = -(ones(ndv,1)*normpdf(-betas')).*dbetas;
                dPt = sum(dPis,2);
                %dbeta = - dinvnc(Pt)*dPt;
                %d(logP) = 1/P*dP = f(b)/P*db; db = dP/(f(b))
                if isinf(beta)
                    dbeta = -dPt;
                else
                    dbeta = -dPt/normpdf(-beta);
                end
            end
        end
    end
    %vargout={Fs,Us,Sigs,en,vol,T,floc,u_flamb,Lambd,D,freq,dF,du,dsig,den,
    %dvol};

else
    %Deterministc case

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nonlinear Analysis
    if HYPLAS_flag
        typeRB = -RBtype;
        if typeRB==0, typeRB = -3; end
        [~,Us,Sigs,en,vol] = Hyplas_analysis( xvalu, command_string, typeRB);

        % Sensibility - Finite Difference Computation
        ndof=length(Us(:,1));
        dx=(norm(xvalu)+1)*1e-5;
        du = zeros(ndvab,ndof);
        dsig = zeros(ndvab,size(Sigs,1));
        den = zeros(ndvab,1);
        dvol = zeros(ndvab,1);
        for i=1:ndvab
            x1=xvalu;
            x1(i) = x1(i) + dx;
            [~,Us1,Sigs1,en1,vol1  ] = Hyplas_analysis( x1, command_string, typeRB);
            du(i,:) = (Us1(:,end)-Us(:,end))'/dx;
            dsig(i,:) = (Sigs1(:,end)-Sigs(:,end))'/dx;
            den(i,1) = (en1-en)/dx;
            dvol(i,1) = (vol1-vol)/dx;
        end
        nmdfl=1;
        lamb=zeros(nmdfl,1);    D=[];    freq=[];floc=[];
        area=0;


    else
        %[Fs,Us,Sigs,en,vol,T]=NLfesol(vargs{:});
        [~,Us,Sigs,en,vol,~,floc,~,lamb,~,~,~,du,dsig,den,dvol]=NLfesol(vargs{:});
        %lamb=zeros(nmdfl,1);    D=[];    freq=[];floc=[];
        dsig=dsig';
    end
    
    u=Us(:,end);
    sig=Sigs(:,end);

end
        
if (Stochast_moment)&&(Stochast_on)
    %[fre,gre,fob,gob]=totim2D(tres,tobj,vol,dvol,sn,dsn,sig,dsig,u,du);
    [fre,gre,fob,gob]=Stoch_TOtim(tres,tobj, v_out0, M_out, S_out, dv_out0, dM_out, dS_out);
%     fob = vol;
%     gob = dvol;
%     fre=[abs(sig) + 3*ssig];
%     gre=[(dsig'.*([1;1;1]*sign(sig')))];% + 3*sdsig')
    if strcmp(Method,'R2DO')||strcmp(Method,'R2DOmc')
        fre=[fre; ReliabC];
        gre=[gre dReliabC];
    end
else
    nel=length(sig);
    nmdfl=1;%number of flamb modes
    dfloc=zeros(ndvab,nel);dlamb=zeros(ndvab,nmdfl);
    [fre,gre,fob,gob]=totim(tres,tobj,vol,dvol,en,den,sig,dsig,u,du,floc,dfloc,ndvab,area,props,lamb(1:nmdfl),dlamb(:,1:nmdfl));
    if Stochast_on
        fre=[fre; ReliabC];
        gre=[gre dReliabC];
    end
end

% NORMALIZE Func Objective
global noinitial_design Fob0 
noinitial_design = 0; % Skip normalization (comment line to use it!)
if isempty(noinitial_design)
    Fob0 = fob;
    noinitial_design = 1;
    fob = fob./Fob0;
    
elseif noinitial_design == 1;
    fob = fob./Fob0;
end

%disp([xvalu, fob])


