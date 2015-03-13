function [G,dG,failure_sets] = Reliab_approx(reliaboptimiz,M_out,dM_out,S_out,dS_out)

global betatarg G_limits set_Failure_FORM

% Setting Reliablt Analss options
%uniq_Failure=0;
 uniq_Failure=1;
%tol_Failure=0.01; 
% tol_Failure=5;



%PMA - Performance mensure approach
%M+b_lim.S < G_limit

%RIA - Reliability index approach
%M+beta.S = G_limit, beta>b_lim
% reliaboptimiz = 'RIA';
% Precomutating!
    ABS_Failure_func = 0;%Regular formulation
%     ABS_Failure_func = 1;%Paper error formulation case
    if (sign(G_limits)==-1)||(ABS_Failure_func)
        g(:,1)=M_out; 
        dg(:,:,1) = dM_out;
    else
        g(:,1)=abs(M_out);
        sgn = sign(M_out(:)');
        oneM=dM_out(:,1)*0+1;
        dg(:,:,1) = (oneM*sgn).*dM_out;
    end
    
    m = g;
    v = S_out.^2;
    mu = log((m.^2)./sqrt(v+m.^2));
    sigma = sqrt(log(v./(m.^2)+1));
    PF=1-cdf('lognormal',abs(G_limits),mu,sigma);
    betap=-norminv(PF);
    if strcmp(reliaboptimiz,'RIA')
        %M+beta.S = G_limit, beta>b_lim
        %beta = (G_limit-M)/S
        %dbeta= -dM/S + M/S.dS/S
        %dbeta= (M.dS-dM.S)/S^2
        G  =  (abs(G_limits) - g)./S_out;
        G  = min(G,betap);
        vabs=ones(size(dM_out,1),1);
        dG = (dg.*(vabs*S_out') - (vabs*G').*dS_out)./(vabs*S_out').^2;
    elseif strcmp(reliaboptimiz,'PMA')
        %M+b_lim.S - G_limit < 0
        g1(:,1)    =  g(:,1)  +betatarg* S_out;
        dg(:,:,1) = dg(:,:,1)+betatarg*dS_out;
        %M+s*beta=Glim, s=(abs(G_limits) - g)./betap;
        %M+s*btar=G2  ,G2=g + btar/betap*(abs(G_limits) - g)
        s=(abs(G_limits) - g)./betap;
        g2 = g + betatarg*s;
        g  = max(g1,g2);
        G = g./abs(G_limits) - sign(G_limits);        
        dG = dg/abs(G_limits);
    end
    

%         % Computing all failure function of the mean vab values
%         set_Failure_FORM=0;
%         [gV,dgU] = feval(G_fun,Mu,varargin{:});

        if uniq_Failure
            % Select worst cases failue functions
            if strcmp(reliaboptimiz,'RIA')
%                 failure_sets=find((G<tol_Failure)&(G>1))';
%                 if find(G<1,1)
%                     failure_sets=[];
%                 end
                [~,failure_sets] = max(G);
            else
                % PMA
%                 disp('need implementation... [Raliab_approx]')
%                 failure_sets=find((G<1)&(G>-1))';
                [~,failure_sets] = max(G);
            end
        elseif uniq_Failure==0
            % Select all failue functions
            n_set_failure=length(gV);
            failure_sets=1:n_set_failure;
            if strcmp(reliaboptimiz,'RIA')
                failure_sets=find((G<tol_Failure)&(G>1))';
                if find(G<1,1)
                    failure_sets=[];
                end
            elseif strcmp(reliaboptimiz,'PMA')
                % PMA
%                 disp('need implementation... [Raliab_approx]')
                failure_sets=find((G<1)&(G>-1))';
            end
        end
        
        if ~isempty(set_Failure_FORM)
            if set_Failure_FORM~=0
                %failure_sets = set_Failure_FORM;
                failure_sets = 1;
                G = G(set_Failure_FORM);
                dG = dG(:,set_Failure_FORM);
            end
        end
%         n_fail_func = length(failure_sets);