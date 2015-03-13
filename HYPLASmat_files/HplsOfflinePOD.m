function [D,F,Cs,Toff,x0s,niters,Dsig,Den,Dvol]=HplsOfflinePOD(n,arq,funame,varargin0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Offline Stage - Basis Samples Generate %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global offline_stage echo Den Dvol Dsig
global Stochast_on RandVarDist Correl 
global convergence_RBM i_convR f_convR savedState
try
    if offline_stage
        % Ensure erro (force offline)
        load xxxOffDate
    else
        load OffDate
    end
    
catch %#ok<CTCH>
%     offline_stage=1; %#ok<NASGU>
    T0=cputime;
    echo0=echo;
    echo=1; %#ok<NASGU>
    %Zdu=[];

	ZKs={};    
        
    % Fixing seed
    defaultStream = RandStream.getDefaultStream;
    defaultStream.State = savedState;
    
    if Stochast_on
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reliability Analysis
        
        ft = RandVarDist(1,:);
        Mu = RandVarDist(2,:);
        Su = RandVarDist(3,:);
        Cov= Correl.*(Su'*Su);
    
        Vdim=size(Mu);
        [Nv,ic]=max(Vdim); %#ok<NASGU>
        Unc=zeros(n,Nv);
        for i=1:Nv
            switch ft(i)
            case 1 %PDF Normal
                Unc(:,i) = lhsnorm(Mu(i), Cov(i,i), n);
            case 2 %PDF LogNormal
                m=Mu(i);
                s=Cov(i,i);
                MU = log(m^2 ./ sqrt(s+m^2));
                SIG = sqrt(log(s/m.^2 + 1));
                Unc(:,i) = lognrnd(MU, SIG, n,1);
            end
        end
        Xrv=Unc;
    else
        Xrv=[];
    end
    
    %% Desig Variables Distribution
    Lx = arq.xmax-arq.xmin;
    nv=length(arq.x0);
    x0s=lhsdesign(n,nv)*Lx+arq.xmin;
    if convergence_RBM
        x0s=x0s(i_convR:f_convR,:);
        Xrv=Xrv(i_convR:f_convR,:);
    end
    
    %% Output Basis Computation
    nR=size(x0s,1);
    D=[];Den=[];Dvol=[];Dsig=[];%nitt=0;
    F=cell(nR,1);niters=zeros(nR,1);
    %area = arq.x0(arq.arex);
    %varargin0={area,props,fext,glb,tipo,-1,0};
    for i=1:nR
        %es0=zeros(nel,1);
        varargin_i=varargin0;
        % Design Variables
        varargin_i{1}=x0s(i,arq.arex);
        if Stochast_on
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input Variables Definitions (U)
            varargin_i{4}=Xrv(i,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        fprintf('Off line snapshats: %s \n',num2str([i, x0s(i,:) Xrv(i,:)]))
        %arq.x0=x0s(i,:);
        %area = arq.x0(arq.arex);
        [Fs,Us,Sigs,ien,ivol]=feval(funame, varargin_i{:});
%         nitt0=sum(niters);
        niters(i)=size(Us,2);
        %for j=1:niters(i), ZKs{nitt0+j}=Ks{j}; end
        %ZKs={ZKs{:},Ks{:}};
        Den=[Den ien];Dvol=[Dvol ivol];
        D=[D Us];
        Dsig=[Dsig Sigs];
        F{i}=max(abs(Fs))';
    end
    
    %% Linear Dependence
    lld=0;
    while ~isempty(lld)
    cit=[0; cumsum(niters)]+1;
    lld=lindep(D');
    if ~isempty(lld)
        disp('Z linearmente dependente, colunas:')
        disp(num2str(lld))
        nlld=size(lld,1);k_lc=zeros(nlld,2);
        for i=1:nlld
            ill=lld(i,2);
            k=find(cit>ill,1);
            dk=lld(i,2)-cit(k-1)+1;
            %kk=cit(k)-dk;
            k_lc(i,:)=[k-1 dk];
            %niters(k-1)=niters(k-1)-1;
        end
        for i=1:nR
            %kl=(k_lc(:,1)==i);
            %F{i}(k_lc(kl,2))=[];
        end
        %[aa,bb]=sort(k_lc(:,1));
        %k_lc(bb,:);
        %F{k-1}(dk)=[];
        %cit(k:end)=cit(k:end)-1;
        D(:,lld(:,2))=[];
        
        %k_lld=lld(:,2);
        %F(:,lld(:,2))=[];
    end
    end
    
% lld=0;
%     while ~isempty(lld)
%     lld=lindep(Zdu');
%     if ~isempty(lld)
%         disp('Zdu linearmente dependente, colunas:')
%         disp(num2str(lld))
%         Zdu(:,lld(:,2))=[];
%     end
%     end

%% K aproximation
% K*Z.a=R
% (K0+S[Ki.b])Za=R
% (K0*Z+S[Ki*Z.b])a=R
% (K0n+S[Kin.b])a=Rn
%  [F0_1, F0_2, F0_...]*b = [F0i]
%  [X1_1, X1_2, X1_...      [X1i]
%  [X2_1, X2_2, X2_...      [X2i]
%  [X._1, X._2, X._...      [X.i]
% Cs*b=[F,X]'
% Kni = K0n+S[Kin.b]

    Cs=[];
%     for ir=1:nR
%         for i=1:niters(ir)
%             Cs = [Cs [F{ir}(i);x0s(ir,:)']];
%         end
%     end

%     nitt0=sum(niters);
%     iglnr=1:size(Us,1);
%     for i=1:nitt0
%         Zkn{i}=D'*ZKs{i}*D;
%         %dKn1=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
%     end
%% Saving Data
    T1=cputime;
    Toff=T1-T0
    
    save('OffDate','x0s','D','ZKs','Cs','F','lld','Toff','niters','n')
%     offline_stage=0;
    echo=echo0;
end
