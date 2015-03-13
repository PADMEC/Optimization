function [PF,beta,dbeta,varargout] = Reliab(Mu,Cov,ft,G_fun,Method,DinamcVar,Trans_type,varargin)
% Reliability analysis
% [PF,beta] = Reliab(Mu,Cov,ft,G_fun,Method,DinamcVar,Trans_type)
% Compute probabilite of failure (PF) and Reliability index (beta)
% of a given failure function (G_fun) and random variable U of mean Mu,
% Covariance matrix (Cov) and pdf (ft). 
% Uses MC or FORM method: Method='MC' or Method='FORM'
% DinamcVar, can be num of integration point or fist guess of MPP
%   
% Example of input data:
%     Mu=[1500,13];
%     Su=[200.,2.];
%     C=[1, .6;.6, 1];
%     C=[1, .0;.0, 1];
%     Cov = C.*(Su'*Su);
%     ft=[2,2];
%     DinamcVar = 30000;Method='MC';
%     DinamcVar = 0;Method='FORM';
%     %DinamcVar = 0;Method='FERUM';
%     G_fun = 'g_constr';
% 

global HYPLAS_flag
if ~HYPLAS_flag
    [Mu,Cov]=check_DesRand(Mu,Cov,varargin{:});
end

global linStressComputed var_type echoReliab
global savedState set_Failure_FORM failure_sets
defaultStream = RandStream.getDefaultStream;
defaultStream.State = savedState;

        
    Su=diag(Cov)'.^.5;
    %Cov = C.*(Su'*Su);
    C=Cov./(Su'*Su);
    
    if nargin<7
        Trans_type=2;
    end

    Vdim=size(Mu);
    [Nv,ic]=max(Vdim);
    
    if linStressComputed==1
        linStressComputed=0;
        var_type = 1;
        [g,dg0]=feval(G_fun,zeros(0,Nv),varargin{:});
        linStressComputed=1;
    end
    var_type = 2;
    
    if strcmp(Method,'MC')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MC
        n=DinamcVar;
        
        [U,V]=Generate_MC_points(Mu, Cov, ft, n,Trans_type);
        
        %U=Unc;
        % Defining
        [~,dgx0] = feval(G_fun,Mu,varargin{:});
        nv=length(dgx0);
        g=zeros(1,n);
        dgx=zeros(nv,n);
        
        % Computing values
        var_type=1;
        for i=1:n
            [g(i),dgx(:,i)] = feval(G_fun,U(i,:),varargin{:});
            if mod(i,500)==0
                fprintf('MC i: %d , mean (G): %d\n',i,mean(g(1:i)));
            end
        end
        var_type=0;
        
        % Probability of failure
        PF = sum(g>0)/n;
        if PF==0
            mg=mean(g);
            sg=std(g,1);
            %m+x*s=0;
            x=-mg/sg;
            PF = 1-normcdf(x);
            PF(PF>1/n)=1/n
        end
        
        beta=-norminv(PF);
        dbeta=mean(dgx');
        if echoReliab, disp([beta, PF]), end
        
        % Verify MC
        if (Nv==2)&&(echoReliab)
            figure, scatter(U(:,1),U(:,2),40,g,'filled')
            hold on, plot(U(g>0,1),U(g>0,2),'.k')
%             figure, plot(U(:,1),U(:,2),'o')
%             hold on, plot(U(g>0,1),U(g>0,2),'.k')
            %figure, scatter3(U(:,1),U(:,2),U(:,3),40,g)
            %hold on, plot3(U(g>0,1),U(g>0,2),U(g>0,3),'.k')
            %y=6.5:.1:29;x=1./(1/5-1./y);
            %hold on,plot(x,y), 
            plot(Mu(1),Mu(2),'+r')
            title([beta, PF])

            %V=(U-ones(n,1)*Mu)*Juv;

            %V=V1;
            %nn=length(x);
            %X=([x',y']-ones(nn,1)*Mu)*Juv;
            figure, scatter(V(:,1),V(:,2),40,g,'filled')
            hold on, plot(V(g>0,1),V(g>0,2),'xk')
            plot(0,0,'+r')
            %hold on,plot(X(:,1),X(:,2)), plot(0,0,'+k')
            title([beta, PF])
        elseif (Nv==3)&&(echoReliab)
            figure, scatter3(U(:,1),U(:,2),U(:,3),40,g,'filled')
            hold on, plot3(U(g>0,1),U(g>0,2),U(g>0,3),'.k')
            title([beta, PF])
            
        end
        
    elseif strcmp(Method,'FORM')||strcmp(Method,'R2DO')
        
        % Rand Vab Transformation - Nataf
        [parmeters,Jzu,Juz,S1] = preproc_randomvar(Mu,Cov, ft,Nv, Trans_type);
        
%         betas = failure_sets*0;
        
        failure_seted = 0;
        if ~isempty(set_Failure_FORM)
            if set_Failure_FORM~=0
                failure_sets = set_Failure_FORM;
                failure_seted = 1;
            end
        end 
        % Run Reliability analisys of the set of failue functions
        k_form=0;
        for i_reliab=failure_sets
            k_form=k_form+1;
            set_Failure_FORM=i_reliab;
            % Run FORM
            FORM_Reliab %(DinamcVar,parmeters,Jzu,Juz,S1,G_fun,Mu,varargin{:})
            betas(k_form)=beta;
            dPis(:,k_form) = -normpdf(-beta)*dbeta;
            dbetas(:,k_form)= dbeta;
        end
%         
%         betas(betas==0)=inf;
%         % Compute failure probabilities
        Ps = normcdf(-betas );
        Ps(isnan(Ps))=0;
        % Compute union of failure probabilities
        PF = sum(Ps);
%         beta = - inv_norm_cdf(PF);
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Sens
%         dPt = sum(dPis,2);
%         %dbeta = - dinvnc(Pt)*dPt;
%         %d(logP) = 1/P*dP = f(b)/P*db; db = dP/(f(b))
%         if isinf(beta)
%             dbeta = -dPt;
%         else
%             dbeta = -dPt/normpdf(-beta);
%         end

        beta=betas;
        dbeta=dbetas;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~failure_seted
            set_Failure_FORM=[];
        end
    
    elseif strcmp(Method,'FERUM')
        % FERUME
        % Marginal distributions for each random variable
        % probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];
        probdata.marg =[ft(:),Mu(:), Su(:), Mu(:)];
        probdata.correlation = C;
        gfundata.expression=G_fun;
        gfundata(1).flag_sens  = 0;
        analysisopt.echo_flag=1;

        %gfundata.evaluator= 'basic';
                  %gfundata.type= 'expression'
        [ formresults, probdata ] = form(1,probdata,analysisopt,gfundata,[],[],varargin{:});

        formresults.beta        
        PF=acu_aprx(-formresults.beta);
        beta = formresults.beta;
        disp([formresults.beta, PF])
    end

end

function df = FiniteDiff(fun,x,F0,varargin)
global grad_calculation
    grad_calculation = 1;
    %grad_calculation = 0;
    nx = size(x);
    nv = nx(2);
    np = nx(1);
    if (nx(2)==1)
        x=x';
        nv = nx(1);
        np = nx(2);
    end
    
    %x = [x1,x2..];
    de = x/1e4;
    de(de==0)=1e-5;
    %F0 = feval(fun,x,varargin{:});
    df = zeros(nv,np);
    xi = zeros(np,nv);
    %for ip=1:np
    ip=1:np; %FOR
        for iv=1:nv
            
            xi(ip,:) = x(ip,:);
            xi(ip,iv) = xi(iv) + de(iv);
            
            F1 = feval(fun,xi,varargin{:});
            
            df(iv,ip) = (F1-F0)/de(iv);
        end
    %end
    
    grad_calculation = 0;
end


function Y=acu_aprx(y)

    p=[+0.23164190
    +0.31938153
    -0.35656378
    +1.78147794
    -1.82125598
    +1.33027443];

    w=1/(1+p(1)*norm(y));
    z=0;
    for i=6:-1:2
        z=w*(z + p(i));
    end
    
    if y>0
        Y = 1 - z.*exp(-y.^2/2)/(2*pi)^.5;
    else
        Y =     z.*exp(-y.^2/2)/(2*pi)^.5;
    end
end



function pointest(U,Mu,Su,Cu,ft,par)
global echoReliab
    lm=abs(Mu);lm(lm==0)=1;
    ls=abs(Su);ls(ls==0)=1;
    me=norm((mean(U)-Mu)./lm);
    ms=norm((std(U)-Su)./ls);
    n=size(U,1);
    Mc=((U-ones(n,1)*Mu)'*(U-ones(n,1)*Mu))/n;
    mc=norm(Mc - Cu)/norm(Mc);

    disp([me, ms, mc])
    nv=size(U,2);
    for i=1:nv
        if echoReliab
            figure,hold on, 
            plot(sort(U(:,i)),(1:n)/n,'.-k')
        end
        switch ft(i)
        case 1  % Normal distribution
            U(:,i) = Z(:,i)*s+m;
            Y=cdf('norm',sort(U(:,i)),Mu(i),Su(i));
            %U(:,i) = Z(:,i)+m;
        case 2  % Lognormal distribution
            Y=cdf('lognorm',sort(U(:,i)),par(i,1),par(i,2));
        end
        if echoReliab
            plot(sort(U(:,i)),Y,'.-r')
        end
    end

end