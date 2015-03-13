function [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc( Mu,Cov,ft,S_fun,Method,N1,varargin)
%[ v_out0, M_out, S_out ] = StochastiCalc(munew, Mu,Cov,ft,fun,Method,DinamcVar,1,vargs{:})
% Comput statistical moment (mean and std) of a function fun
% fun: [v_out]=Stoch_Analysis(x,varargin)

%global N1
%global mont_carl W_PCM

global var_type test_importance sens ihistogram HYPLAS_flag
% Derivative with respect to the Design Variables
var_type = 1;

% Derivative with respect to the Ramdom Variables
if test_importance
    var_type = 2;
end


% Scatter PLOTs! 2D or 3D
LOGC_PLOTEST=0;
% HISTOGRAM PLOT!
% ihistogram=1;
% if ihistogram||LOGC_PLOTEST
%     sens=0;
% else
%     sens=1;
% end

if strcmp(Method,'R2DOmc')
    %PCM, n = 2
    N1=2;
end

if ~HYPLAS_flag
    [Mu,Cov]=check_DesRand(Mu,Cov,varargin{:});
end

% NonLinear Analysis outputs
%vargout={Fs,Us,Sigs,en,vol,T,floc,u_flamb,Lambd,D,freq,dF,du,dsig,den,dvol

[vargout]=feval(S_fun,Mu,varargin{:});
%vargout=Stoch_Analysis();

%% Multiple output computation
if iscell(vargout)
    n_out=length(vargout);
    if sens
        n_grad = n_out/2;
        n_out = n_grad;
        dv_out0={vargout{n_out+(1:n_grad)}};
        nvar=size(dv_out0{1},1);
    else
        n_grad = 0;
        dv_out0=0;
    end
    v_out0={vargout{1:n_out}};
else
    n_out=1;
    v_out0={vargout};
    n_grad=0;
    dv_out0=[];
    %nvar=0;
end
% n_out=4;
% v_out = {vol,u,vonmis,sn};

%% Defining Outputs
M_out=cell(1,n_out);
F2_out=cell(1,n_out);
S_out=cell(1,n_out);
dM_out=cell(1,n_grad);
dS_out=cell(1,n_grad);
FdF_out=cell(1,n_grad);

for i_out = 1:n_out
    M_out{i_out}= zeros(size(v_out0{i_out}));
    F2_out{i_out}=M_out{i_out};
end

for i_out = 1:n_grad
    dM_out{i_out}= zeros(size(dv_out0{i_out}));
    FdF_out{i_out}=dM_out{i_out};
end

%% Generating Integration Points
if strcmp(Method,'MC')
    % MC
    Trans_type=1;
    [U]=Generate_MC_points(Mu, Cov, ft, N1,Trans_type);
else
    % PC
    [U,W_PCM]=Generate_PC_points(Mu, Cov, ft, N1);
end

if LOGC_PLOTEST==1
 figure, plot(U(:,1),U(:,2),'.')
    nrv = length(Mu);
    if strcmp(Method,'PC')
        figure
        if nrv>=2
            scatter3(U(:,1),U(:,2),prod(W_PCM,2),50,prod(W_PCM,2),'filled')

        else
            plot(U,W_PCM,'.')
        end
    end
%     
%     for i_out = 1:n_out
%         Hfig(i_out)=figure;title(num2str(i_out));hold on
%     end
%     HSN=figure;hold on,
%     HUi=figure;hold on,
%     HMaxu=figure;hold on,
%     HSi=figure;hold on,
%     HMaxS=figure;hold on,
    
end

np=size(U,1);
if ihistogram==1
    MC_Data=cell(1,n_out);
    %MC_Data=zeros(n_out,np);
end

%% Begining Point Evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MC/PC CALCULATIOS %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for imc=1:np
        
        %DESIGN VARIABLES
        x=U(imc,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Point Analysis
        [vargout]=feval(S_fun,x,varargin{:});
        % for multiples output function

        if iscell(vargout)
            v_out={vargout{1:n_out}};
            dv_out={vargout{n_out+(1:n_grad)}};
        else
            n_out=1;
            v_out={vargout};
            n_grad=0;
            dv_out=[];
            nvar=0;
        end
        
        %$ STACHASTICTICS $%
        if strcmp(Method,'MC')
            % MC
            fprintf('i = %d, x = %s \n',imc, num2str(x))
            W_m = 1/N1;
            W_s = 1;
            
        elseif strcmp(Method,'PC')||strcmp(Method,'R2DO')||strcmp(Method,'R2DOmc')
            % PC
            W_m = W_PCM(imc);
            W_s = W_PCM(imc)*N1;
        end
        
%         for i_out = 1:n_out
%             %M(f)=SUM[fi]
%             M_out{i_out}= M_out{i_out} + v_out{i_out}/N1;
% 
%             %S(f)=SUM[fi^2]
%             F2_out{i_out} = F2_out{i_out} + v_out{i_out}.^2;
%         end

        for i_out = 1:n_out
            %PC_mean=w*Zpc
            M_out{i_out}= M_out{i_out} + v_out{i_out}*W_m;

            %PC_std=(W*Zpc.^2-pc_mean^2)^.5
            F2_out{i_out} = F2_out{i_out} + v_out{i_out}.^2*W_s;
        end

        if sens
            for i_out = 1:n_grad
                %COMPUT dM_out = dFi and FdF_out = Fi*dFi ...
                dM_out{i_out}= dM_out{i_out} + dv_out{i_out}*W_m;
                %for ivar=1:nvar
                    dFdFi = (ones(nvar,1)*v_out{i_out}'*W_s).*dv_out{i_out};
                    FdF_out{i_out}= FdF_out{i_out} + dFdFi;
                %end
            end
        end
        % HISTOGRAM PLOT
        if (ihistogram==1) || (LOGC_PLOTEST==1)
            for i_out = 1:n_out
                %PC_mean=w*Zpc
                MC_Data{i_out}(:,imc)=v_out{i_out};
            end
        end
        
    end % End evaluations
        
        
    if LOGC_PLOTEST==1
        try
            % displaciment
            MCdatas(:,1) = MC_Data{2}(661,:);
            % Max mean sig
            [maxval,im]=max(mean(abs(MC_Data{3}')));
            MCdatas(:,2) = MC_Data{3}(im,:);
            % Max std sig
            [~,imsd]=max(std(abs(MC_Data{3}')));
            MCdatas(:,3) = MC_Data{3}(imsd,:);
            % Saving
            arqname=['Statistc_' num2str(np) '_x1_' num2str(varargin{1}(1))];
            save(arqname,'MCdatas','U')
            disp('OK!! saved')
            arqname
        catch
            save(['Statistc_' num2str(np) '_x1_' num2str(varargin{1}(1))],'MC_Data','U')
            disp('not OK!! saved')
        end
        for i_out = 2:3
%                 Hfig(i_out)=figure;hold on
            im=1;
            if (nrv<=3)&&(nrv>=2)
            if size(MC_Data{i_out},1)>1
%                 [maxval,im]=max(abs(MC_Data{i_out}));
                [~,im]=max(mean(abs(MC_Data{i_out}')));
                
                [~,imsd]=max(std(abs(MC_Data{i_out}')));
                maxval=MC_Data{i_out}(imsd,:);
                figure
                scatter3(U(:,1),U(:,2),...
                 U(:,nrv),50,maxval,'filled')
                title(['max SD ' num2str(i_out)])
%                 figure,plot(im,'.')
            end
            end
            maxval=abs(MC_Data{i_out}(im,:));
            
            figure,hist(maxval,30);
            title(num2str(i_out))
            
            %set(0,'CurrentFigure',Hfig(i_out))
            if (nrv<=3)&&(nrv>=2)
                figure,hold on,title(num2str(i_out))
                grid on,xlabel('U1'),ylabel('U2')
                if i_out==2
                    title('d')
                elseif i_out==3
                    title ('Max (\sigma)')
                end
            end
            if nrv==2
                scatter3(U(:,1),U(:,nrv),...
                 maxval,50,maxval,'filled')
                zlabel('F_out')
            elseif nrv==3
                scatter3(U(:,1),U(:,2),...
                 U(:,nrv),50,maxval,'filled')
                zlabel('U3')
            end
            colorbar
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X0=[xsmc; randcoordsmc; randpropsmc; randbondsmc]
 %   [p,xi,I,regpart,fi]=aproximateqncs(xsmc,dp,xs,g,f);
%    Q_MC()
%     [xi,regpart,fi]=Q_MC(xsmc,xsnmc,munew,sdevi);
% 

if N1>30
    RN=1/(N1-1);
else
    RN=1/(N1);
end

for i_out = 1:n_out
    %PC_std=(W*Zpc.^2-pc_mean^2)^.5
    S_out{i_out} = ((F2_out{i_out}-N1*M_out{i_out}.^2)*RN).^.5;
end


if sens
    for i_out = 1:n_grad
        %if S_out{i_out}.^2)<1e-9, dS_out{i_out}=zeros(nvar,1); else

            %d(std)/dx = 1*/S*(fi*W*dfi + NMdM) = 
            Si = ones(nvar,1)*S_out{i_out}';
            no_S = (S_out{i_out}.^2)<1e-12;
            Si(:,no_S)=1;

            % S = (Fi.^2-Fm^2)^.5
            % ds = .5(Fi.^2-Fm^2)^-.5*(2*FidFi-2FmdFm)
            % ds = (FidFi-FmdFm)/S
            NMdM = N1*(ones(nvar,1)*M_out{i_out}').*dM_out{i_out};
            dS_out{i_out}=(FdF_out{i_out}-NMdM)*RN./Si;
        %end
    end
end
% TEST APROXIMAÇÃO %
%[Amean, Astd]=Surfsthocalc(xsmcp,pmcp,snmc)
%%[muhat,sigmahat] = normfit(snmc);
%[Amean, Astd]=my_sthocalc(snmc,pmcp)
%[Amean, Astd]=sthocalcMSdt(xsmcp,pmcp,snmc)
%%%%%%%%%%%%%%%%%%%%
ihistogram=1;
if ihistogram==1
    for i_out = 1:n_out
%         [maxval,im]=max(abs(M_out{i_out}));
%         figure
%         hist(maxval,40)
%         title(num2str(i_out))
    end
    %figure, hist(MC_Data{2}(8,:),50)
    %figure, hist(MC_Data{3}(264,:),50)
%     figure
%     plot(snmc,pmcp,'.')
%     figure, plot(snmc,'.'),hold on
%     plot([1 length(snmc)],[sn0 sn0],'g')
%     plot([1 length(snmc)],[Msn Msn],'k')
%     plot([1 length(snmc)],[Msn+Ssn Msn+Ssn],'r')
%     plot([1 length(snmc)],[Msn-Ssn Msn-Ssn],'r')
end