function [Uopts,Gopts,dgUs,Converg]=PMA_RBDO(Mu,Cov,ft,G_fun,varargin)
% Reliability constraint

global betatarg var_type failure_sets set_Failure_FORM echoReliab HYPLAS_flag
    var_type = 2; % Random Variables Space - [Optimizations
    Trans_type = 2;

    if ~HYPLAS_flag    
        [Mu,Cov]=check_DesRand(Mu,Cov,varargin{:});
    end
    nrv=length(Mu);
    [parmeters,Jzu,Juz,S1] = preproc_randomvar(Mu,Cov, ft,nrv,Trans_type);
    Su=diag(Cov)'.^.5;
%     [C,dC]=PMA_reliab_constr(V,Mu,Su,Jzu,ft,nrv,parmeters,G_fun,varargin{:});
%     [gV,dgV]=PMA_reliab_obj(V,Mu,Su,Jzu,ft,nrv,parmeters,G_fun,varargin{:});
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
failure_seted = 0;
if ~isempty(set_Failure_FORM)
    if set_Failure_FORM~=0
        failure_sets = set_Failure_FORM;
        failure_seted = 1;
    end
end 

i_set=0;
% dgUs=zeros(nrv,length(failure_sets));
Uopts=zeros(length(failure_sets),nrv);
Gopts=zeros(size(failure_sets));
Converg = zeros(1,nrv);

maxiter=20;
for i_fail=failure_sets;
    set_Failure_FORM = i_fail;
    i_set=i_set+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ITERACTIVE METHOD!!!!!!!!!!!!!
    U = Mu;V = zeros(nrv,1);
    iter=0;notconverg=1;g=-inf;
    while notconverg
        iter=iter+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computing V initial
        % V = dgV*beta/norm(dgV)
        % U = S.V + M
        g0=g;
        [g,dgU] = feval(G_fun,U,varargin{:});
        Jzv = inv_normeq(Mu,Su,U,ft);
        dgV = Jzu*Jzv*dgU;
        
        unidv = dgV/norm(dgV);
        % Convegence
        rest=V'*unidv-betatarg;
        
        % New point
        V = unidv*betatarg;
        U = transNf(V,Jzu,ft,nrv,parmeters);
        
        % Check Convegence
        notconverg=(abs(rest)>1e-7)&&(iter<=maxiter);
        if echoReliab 
            if g<g0
                fprintf('Iter %d, G0: %g, G1= %g \n', iter, g0, g)
            end
            fprintf('Iter %d, G Constraint: %g, Convergence = %g \n', iter, g, rest)
            if ~notconverg
                disp(U)
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post-proc
    var_type = 1; % Design Variables Space - Optimizations
    [Gopt,dgU] = feval(G_fun,U,varargin{:});
    var_type = 2; % Design Variables Space - Optimizations
    Uopts(i_set,:)=U;
    Gopts(i_set)=Gopt;
    dgUs(:,i_set)=dgU';
    Converg(i_set)=iter>=maxiter;
    

%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % TEST Obj Gradient
%             [gV0,dgV0]=PMA_reliab_obj(V0,Mu,Su,Jzu,ft,nrv,parmeters,G_fun,varargin{:});
%             U0 = transNf(V0,Jzu,ft,nrv,parmeters);
%             [gU0,dgU0] = feval(G_fun,U0,varargin{:});
%             for i=1:nrv
%                 dx = zeros(nrv,1);
%                 dx(i)=V0(i)/10000;
%                 [gV1,dgV1]=PMA_reliab_obj(V0+dx,Mu,Su,Jzu,ft,nrv,parmeters,G_fun,varargin{:});
%                 Grads(i)=(gV1-gV0)/dx(i);
%                 
%                 dx = zeros(nrv,1);
%                 dx(i)=U0(i)/10000;
%                 [gU1,dgU1] = feval(G_fun,U0+dx',varargin{:});
%                 GradU(i)=(gU1-gU0)/dx(i);
%             end
%             [Grads;dgV0']
%             [GradU;dgU0']
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Mathematical programing optimization!!
%     options=optimset;
%     %options=optimset(options,'LargeScale','on');'SubproblemAlgorithm','cg'
%     options=optimset(options,'GradObj','on','GradConstr','on','MaxFunEvals',40000,'Algorithm' ,'sqp');%'active-set');%'interior-point');%,
%     options=optimset(options,'MaxIter',100,'Display','iter','DerivativeCheck','off', 'LargeScale','off');%,'MaxTime',9);
%     %options=optimset(options,'LargeScale','off','Algorithm','interior-point');
%     vb=[-betatarg;betatarg]*ones(1,nrv);
%     
%     [Vopt,Gopt,Converg,Resumo]=fmincon('PMA_reliab_obj',V0,...
%         [],[],[],[],vb(1,:),vb(2,:),'PMA_reliab_constr',options,...
%         Mu,Su,Jzu,ft,nrv,parmeters,G_fun,varargin{:});
%     disp([Gopt-g]/norm(g))
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
if ~failure_seted
    set_Failure_FORM=[];
end
end