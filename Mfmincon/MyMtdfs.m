
%pares=[1 1 0;1 0 1;0 1 1];
pares=gridsamp([zeros(1,nfun);ones(1,nfun)],2);
pares(sum(pares,2)<2,:)=[];
%ones(3)-eye(3);
[probdims,indxprob]=sort(sum(pares,2));
pares=pares(indxprob,:);

% nfsub=sum(pares,2)';
% [nfsub,ksort]=sort(nfsub);
% pares=pares(ksort,:);
% pares=pares(1:3,:);
%npar=find(probdims==3,1)-1;
%nit=(nfun^2-nfun)/2;

%OK
% figure
% grid on
% hold on
% plot3(ft(:,1),ft(:,2),ft(:,3),'.r')
% 
% rand('seed',2481632)
% randn('seed',111111)
logplot=1;

paret_front_dim=2;
b2=betaweigs(paret_front_dim,ptsp,2);

%fis=zeros(size(b2),npar);
f=ft;
diffx=[];
aT=ft;
[TT,mfT,lT]=MO_precomputs(ft,method);

Fpd{1}=ft;
Xpd{1}=x;
pardim=find(probdims==paret_front_dim)';
for ip=pardim
    
    par=find(pares(ip,:));
    inv=find(pares(ip,:)==0);
    %par=kpar(ip,:);
    nf=length(par);
    ff=ft(par,par);
    [T,mf,l,n,t]=MO_precomputs(ff,method);
        
    if probdims(ip)==2
        x0=x(par(1),:);
        if method==5
            x0=x(par(2),:);%comeco das curvas
        end
            x0=x(par(2),:);%comeco das curvas
        tic
        normf=norm(ft(par(2),:));
        if normf<1e-12, normf=1;end
        if norm(ft(par(2),:)-ft(par(1),:))/normf<1e-6
            for i=1:size(b2,1)
                %Point Pareto Curve
				f0(i,:)=ft(par(1),:);
				x(end+1,:)=x0;
				f(end+1,:)=f0(i,:);
                fcount=[fcount 0];
                Cnvrg=[Cnvrg 1];
                Xs(i,:,ip)=x0;
            end
        else
            
%% 2D Pareto Front XXXXXXXXXXX
            aT=ft;
            for i=1:size(b2,1)%:-1:1
                 a=T*b2(i,:)';
                 
                x0=x(end,:);
                % Choose nerest solution
                if (method==6)
                    ax=a.*l'+mf';
                elseif method==5
                    ax=a+mf';
                else
                    ax=a;
                end
                Axs(i,:)=ax;
                [mnd,ix]=min(sum((ax*ones(1,size(f,1))-f(:,par)').^2));
                if norm(ax'-f(end,par))>1.2*norm(ax'-f(ix,par))
                    x0=x(ix,:);
                end

                 
%                 atx=lT+mfT;atx(par)=a';
%                 if method==6
%                     atx=atx.*lT+mfT;
%                 elseif method==5
%                     atx=atx+mfT;
%                 end
%                 aT(end+1,:)=atx;
                %a=a-mf';
                [xt,ffobj,Converg,Resumo]=fmincon(funmtd,[x0 t],[],[],[AVeq t*0],BVeq,...
                    vlb,vub,conmtd,options,fun,con,-1,mf,l,a,n',par,varargin{:});

                funcCount=Resumo.funcCount;
                fcount=[fcount funcCount];
                Cnvrg=[Cnvrg Converg];
                
                
                if method==3||(method==5)
                    t=xt(end);
                    ts(i)=t;
                    xt(end)=[];
                end
                diffx(end+1)=norm(xt-x0);
                x(end+1,:)=xt;
                %x0=xt;
                f0(i,:)=feval(fun,x(end,:),varargin{:});
                
                 options=optimset(options,'TolCon',1e-4);
                [xi,ffobj,Convrg,Resm]=fmincon('constfun',xt,[],[],AVeq,BVeq,...
                    vlb,vub,'constcon',options,fun,con,inv,f0(i,par),0,f0(i,:),par,varargin{:});
                 options=optimset(options,'TolCon',1e-7);
                if norm(xi-xt)/norm(x0)>1e-4
                    %x0
                    %xi
                    f0(i,:)=feval(fun,xi,varargin{:});
                end
                x(end,:)=xi;
                x0=xi;
                f(end+1,:)=f0(i,:);
                Xs(i,:,ip)=xi;
                
                fprintf('F_%d =   %s\n',size(f,1),num2str(f(end,:),4));
                fprintf('Convergencia: %d   funcCount: %d \n\n'...
                     ,Converg,funcCount);
                %TEST=[-ts'*n+Axs f0]
%                 figure, plot(Axs(:,1),Axs(:,2),'.'), grid on
%                 hold on, plot(f0(:,par(1)),f0(:,par(2)),'.r'), grid on
            end
%% Pos-Proc Pareto Curve
        end
        NBI=toc
        fis(:,:,ip)=f0;
        %OK
        %plot3(f0(:,1),f0(:,2),f0(:,3),'.-k')
%         f=[f; f0];
    end
    %figure
    %scatter3(fi(:,1),fi(:,2),fi(:,3),40,fi(:,3),'filled')
end


    %for i=1:3
    %    maxf(i+1,:)=max(f(:,:,i));
    %end
    %maxf=max(maxf);
    %l=maxf-mf;
    %l=mxf-mf;
    %mnf=l;
%if nfun==4
    ifcoln=nfun;
    if logplot==10
        figure, hold on, 
        maxi=pardim(end);
        for i=pardim%OK
            colr=[ i/maxi (maxi-i+1)/maxi 2*i/maxi-1];colr(colr<0)=0;
           scatter3(fis(:,1,i),fis(:,2,i),fis(:,3,i),40,colr,'filled')
        end
        plot3(ft(:,1),ft(:,2),ft(:,3),'*')
        grid on, colorbar
    end

%%
Fpd{2}=fis;
Xpd{2}=Xs

global mP F_iters Ps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if logplot==1
        if nfun>2
            hplotF=figure,hold on;
            ierotim=find(Cnvrg<1);
            iok=find(Cnvrg>0);
            fok=[f(iok,:)];
            ferr=f(ierotim,:);
                %chin=fill3(ft(:,1),ft(:,2),ft(:,3),'b');
                %alpha(chin,.3);
                %plot3(bt(:,1),bt(:,2),bt(:,3),'.')
                scatter3(fok(:,1),fok(:,2),fok(:,3),50,fok(:,nfun),'filled')
                plot3(ferr(:,1),ferr(:,2),ferr(:,3),'x')
                scatter3(ferr(:,1),ferr(:,2),ferr(:,3),40,ferr(:,nfun))
                colorbar
                grid on, title(['F - ' MOptype])
        end
        end
ind_x=[1:nfun];
%% Pareto Dimension Interations 
while any(probdims==paret_front_dim+1)
    paret_front_dim=paret_front_dim+1

    pre_dim=pardim;
    pardim=find(probdims==paret_front_dim)';

    b=betaweigs(paret_front_dim,ptsp,2);
    aTa=aT;
    
%% Pareto front Boundaries Interations 
    ip=0;
    for ipar=pardim
        %aT=aTa;
        ip=ip+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% DATA %%%%%%%%%%%%%
%% Updating DATA 
        
        par=find(pares(ipar,:));
        inv=find(pares(ipar,:)==0);
        ff=ft(par,par);
        [T,mf,l,n,t]=MO_precomputs(ff,method);

        %%bound solutions index
        surfbond=find(~any(pares(pre_dim,inv),2))';
        
%% Projeting Pareto Boundary    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Incluir todos os otimos %%%%%%%%%% 
%%%%%%%%%% de dimensoes inferiores %%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fs{1}=Fpd{1}(par,par);
        Xpds=Xpd{1}(par,:);
        for i_dim=2:paret_front_dim-1
    
            %%bound solutions index
            prob_dimsi=find(probdims==i_dim)';
            surfbondi=find(~any(pares(prob_dimsi,inv),2))';
    
            Fs{i_dim}=Fpd{i_dim}(:,par,surfbondi);
            
            nbond=length(surfbondi);
            for ibond=1:nbond
                Xpds=[Xpds; Xpd{i_dim}(:,:,ibond)];
            end
            %Fpd{paret_front_dim-1}(:,par,surfbond)
        end
        [Pt,n1,ff,pf,tp,Vr,xvr,P0]=projctpareto(ff,Fs,mf,l);
        
% %% TEST CVT (Plot 1 - boundary)
%         %Space - Vr plan
        if logplot==10
        if (paret_front_dim==3)
             figure, hold on,title('xs Projetct Parameters')
             plot(xvr(1,:),xvr(2,:),'.')
             figure, hold on,title('Projected Points')
             plot3(pf(:,1),pf(:,2),pf(:,3),'.')
%             figure, hold on,
%             plot(xvr(1,:,1),xvr(2,:,1),'.')
%             plot(xvr(1,:,2),xvr(2,:,2),'.k')
%             plot(xvr(1,:,3),xvr(2,:,3),'.r')
%             figure, hold on,
%             plot3(pf(:,1,1),pf(:,2,1),pf(:,3,1),'.')
%             plot3(pf(:,1,2),pf(:,2,2),pf(:,3,2),'.k')
%             plot3(pf(:,1,3),pf(:,2,3),pf(:,3,3),'.r')
        elseif (paret_front_dim==4)
            if ip==1
            figure, hold on,
            plot3(xvr(1,:),xvr(2,:),xvr(3,:),'.m')
            end
%             plot3(xvr(1,:,1),xvr(2,:,1),xvr(3,:,1),'.g')
%             plot3(xvr(1,:,2),xvr(2,:,2),xvr(3,:,2),'.k')
%             plot3(xvr(1,:,3),xvr(2,:,3),xvr(3,:,3),'.r')
%             plot3(xvr(1,:,4),xvr(2,:,4),xvr(3,:,4),'.m')
%             figure, hold on,
%             scatter3(pf(:,1,1),pf(:,2,1),pf(:,3,1),40,pf(:,4,1),'filled')
%             scatter3(pf(:,1,2),pf(:,2,2),pf(:,3,2),40,pf(:,4,1),'filled')
%             scatter3(pf(:,1,3),pf(:,2,3),pf(:,3,3),40,pf(:,4,1),'filled')
        end
        end
%         Space - Vr plan

%% CVT Distributing Points

        npip=size(b,1);
        [xvi0]=generate_intP(xvr',npip);
        
        %X0=b*ff;
        %X0=Pi
        %[X0p,xvi0,tp]=PointProj(Vr,P0,X0);
        if logplot==10
        figure, hold on,title('xs Projetct Parameters')
        if (paret_front_dim==3)
            plot(xvi0(1,:),xvi0(2,:),'*')
            plot(xvr(1,:),xvr(2,:),'.')
        elseif (paret_front_dim>3)
            scatter3(xvi0(1,:),xvi0(2,:),xvi0(3,:),35,...
                xvi0(paret_front_dim-1,:),'filled')
            scatter3(xvr(1,:),xvr(2,:),xvr(3,:),35,...
                xvr(paret_front_dim-1,:),'filled')
        end
        end
        
        [nx0,ndim]=size(xvi0);[ndimf]=size(xvr,2);
        sample_num_cvt=nx0*(ndim+ndimf)*20;
        
        [X,cvt_time]=uniformCVT_IntPoint(xvi0,sample_num_cvt,Vr,xvr,P0);
        %fprintf('CVT - %d Points distributing, time: %f s',nx0,cvt_time)
        
% %% TEST CVT (Plot 2 - Int_point)
%          plotnum(X(:,1),X(:,2))
%          plot(xvi0(1,:),xvi0(2,:),'og')
    if logplot==1
    if (paret_front_dim==4)
        if ip==1
            figure,hold on,
            plot3(xvr(1,:),xvr(2,:),xvr(3,:),'.')  
            plot3(X(:,1),X(:,2),X(:,3),'*')
            plot3(xvi0(1,:),xvi0(2,:),xvi0(3,:),'og')
        end
    elseif (paret_front_dim==3)
        figure,hold on,
        plot(xvr(1,:),xvr(2,:),'.')  
        plot(X(:,1),X(:,2),'*k')
        plot(xvi0(1,:),xvi0(2,:),'.g')
    end
    end

%% Int_point Post-proc

        %auxiliar 'scalar-matrix' vetor
        one_nIP=ones(1,size(X,1));
        
        %Desparametrizing int_points
        p0=X*Vr' + one_nIP'*P0';
        %Desnormalizing int_points
        p0=p0.*(one_nIP'*l)+one_nIP'*mf;

    %Desnormalizing points
        %projected solucions points matrix (aT) and respect x
        %aT=Fpd{1}(par,:);
        %Xa=Xpd{1}(par,:);
        mnPp=size(pf,1);
        MlT=ones(mnPp,1)*l;
        MmfT=ones(mnPp,1)*mf;
        npresub=size(Fpd{paret_front_dim-1},1);
        %atx=zeros(npresub,nfun);ipbond=0;
        %for ibiP=surfbond
        %    ipbond=ipbond+1;
%             if method==6
%                 atx=pf(:,:,ibiP).*MlT+MmfT;
%             elseif method==5
%                 atx=pf(:,:,ibiP)+MmfT;
%             end
            %ind_x=[ind_x nfun+(1:npresub)+(ibiP-1)*npresub];
            %atx(:,par)=pf(:,:,ipbond).*MlT+MmfT;
            %aT=[aT; atx];
            aT=zeros(mnPp,nfun);
            aT(:,par)=pf.*MlT+MmfT;
            Xa=Xpds;

%         figure, hold on,
%             Vplot=ft(par,par);
%             plot3(Vplot(:,1),Vplot(:,2),Vplot(:,3),'.k')
%             Vplot=aT(:,par);
%             plot3(Vplot(:,1),Vplot(:,2),Vplot(:,3),'.')
%         for i=surfbond
%             Vplot=fis(:,par,i);
%             plot3(Vplot(:,1),Vplot(:,2),Vplot(:,3),'.k')
%             %Vplot=pf(:,:,i);
%             %plot3(Vplot(:,1),Vplot(:,2),Vplot(:,3),'.')
%         end
        %Xa=x(ind_x,:);
        
%         for ibiP=pardim
%             atx=fis(:,:,ibiP);
%             aT=[aTa; atx];
%         end
        
%% Plot scripts before sub optimization
        if logplot==1
        if ~isempty(mP)
            %Degenerate Contourn create new point plan
            figure, plot3(mP(:,1)*l(1)+mf(1),...
                mP(:,2)*l(2)+mf(2),mP(:,3)*l(3)+mf(3),'xk')
            scatter3(f(:,1),f(:,2),f(:,3),40,f(:,nfun))
            title('Points to the new project plan')
        end
        end

        for i=pardim%OK
    %       plot3(pf(:,1,i)*l(1)+mf(1),pf(:,2,i)*l(2)+mf(2),pf(:,3,i)*l(3)+mf(3),'xr')
        end

    %     plot3(pf(:,1,1)*l(1)+mf(1),pf(:,2,1)*l(2)+mf(2),pf(:,3,1)*l(3)+mf(3),'xr')
    %     plot3(pf(:,1,2)*l(1)+mf(1),pf(:,2,2)*l(2)+mf(2),pf(:,3,2)*l(3)+mf(3),'xr')
    %     plot3(pf(:,1,3)*l(1)+mf(1),pf(:,2,3)*l(2)+mf(2),pf(:,3,3)*l(3)+mf(3),'xr')

    %     for i=1:size(b,1)
    % 
    %         %ai(:,i)=distpont(b(i,:),b2,f,nfun,apar);
    %         %Vt=(ft(kpar(:,2),:))-(ft(kpar(:,1),:));
    %         %ai(:,i)=distpont2(b(i,:),b2,f,nfun,apar,kpar,Vt,ft);
    %         [a,ta]=Compute_initpoint(b(i,:),b2,pf,nfun,pares,ff,l,mf,tp);
    %         as(:,i)=a-(Pt.*l')*ta;
    %         ai(:,i)=a;
    %         %plot3(ai(1,i),ai(2,i),ai(3,i),'*g')
    %     end
        %p0=ai';
        %figure,plot3(p0(:,1),p0(:,2),p0(:,3),'*g'),grid on
        %OK
        % figure
        % grid on
        % hold on
        % for i=pardim
        %     plot3(fis(:,1,i),fis(:,2,i),fis(:,3,i),'.-k')
        %     %fis=[fis; fis(:,:,i)];
        %     plot3(pf(:,1,i)*l(1)+mf(1),pf(:,2,i)*l(2)+mf(2),pf(:,3,i)*l(3)+mf(3),'xr')
        % end
        % plot3(p0(:,1),p0(:,2),p0(:,3),'x')
        % plot3(ft(:,1),ft(:,2),ft(:,3),'.r')
        % %fis=[ft;fis];
        % plot3(ai(1,:),ai(2,:),ai(3,:),'.')
        %     quiver3(ai(1,:),ai(2,:),ai(3,:),n1(1),n1(2),n1(3),0,'filled');
        %     
        % hold off
        ax=zeros(1,nfun);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Sub optimization
        for i=1:size(b,1)

            a=p0(i,:)';
            %a(par)=p0(i,:)';
            %a(inv)=0;
            %[mnd,ix]=min(sum((a*ones(1,size(f,1))-f').^2));
            %x0=x(ix,:);
            
            %Find int point 
            ax(par)=a;
            [mnd,ix]=min(sum((a*ones(1,size(aT,1))-f(:,par)').^2));
            x0=x(end,:);
            if norm((a'-f(end,par))./l)>2*norm((a'-f(ix,par))./l)
                x0=x(ix,:);
            end 

            disp(i)
            %disp(a')
            if method==4||(method==6)
                a=(a-mf')./l';
            else
                a=a-mf';
            end
            %n=-a'/norm(a);
            %[xt,ffobj,Converg,Resumo]=fmincon('nbifun',[x0
            %t],[],[],AVeq,BVeq,vlb,vub,'nbicon',options,fun,con,-1,mf,l,a,n',par);
            
            [xt,ffobj,Converg,Resumo]=fmincon(funmtd,[x0 t],[],[],[AVeq t*0],BVeq,...
                vlb,vub,conmtd,options,fun,con,-1,mf,l,a,n',par,varargin{:});


            %PLOT NNC CONSTRAINTs
             %hold on, plot3(Ps(:,1),Ps(:,2),Ps(:,3),'.r')
             %hold on, plot3(F_iters(:,1),F_iters(:,2),F_iters(:,3),'.-g')
        %      nn=null(n)*norm(l)/10;
        %      v1=[Ps(1,1) Ps(1,1)+nn(1,1); Ps(1,1) Ps(1,1)+nn(1,2)];
        %      v2=[Ps(1,2) Ps(1,2)+nn(2,1); Ps(1,2) Ps(1,2)+nn(2,2)];
        %      v3=[Ps(1,3) Ps(1,3)+nn(3,1); Ps(1,3) Ps(1,3)+nn(3,2)];
        %      if i==11
        %      Ps=a'.*l+mf;
        %      vv=[Ps; Ps+nn(:,1)'; Ps; Ps+nn(:,2)'];
        %      hold on, plot3(Ps(:,1),Ps(:,2),Ps(:,3),'.r')
        %      hold on, plot3(vv(:,1),vv(:,2),vv(:,3),'r')    
        %      hold on, plot3(F_iters(:,1),F_iters(:,2),F_iters(:,3),'.-g')
        %      FF=F_iters;
        %      PP=Ps;
        %     end
             F_iters=[];
             Ps=[];

            if method==3||(method==5)
                t=xt(end);
                ts(end+1)=t;
                xt(end)=[];
            end
            fcurnt=feval(fun,xt,varargin{:});
            if ~isempty(inv)
%                 options=optimset(options,'TolCon',1e-4);
                [xt,ffobj,Conver,Resum]=fmincon('constfun',xt,[],[],AVeq,BVeq,...
                    vlb,vub,'constcon',options,fun,con,inv,fcurnt(:,par),0,fcurnt,par,varargin{:});
%                 options=optimset(options,'TolCon',1e-6);
            end
            fcurnt=feval(fun,xt,varargin{:});
            
            diffx(end+1)=norm(xt-x0);

            funcCount=Resumo.funcCount;
            fcount=[fcount funcCount];
            Cnvrg=[Cnvrg Converg];
            
            

            x(end+1,:)=xt;
            f(end+1,:)=fcurnt;
            fprintf('F_%d =   %s\n',size(f,1),num2str(f(end,:),4));
            fprintf('Convergencia: %d   funcCount: %d \n\n'...
                 ,Converg,funcCount);
            aT(end+1,:)=ax;
            Xa(end+1,:)=xt;
            disp(fcurnt)
            %disp([x0 ;xt])
            x0=x(end,:);
            xbs(i,:)=x0;
            %fi(i,:)=feval(fun,x(end,:));
            
            fbs(i,:)=fcurnt;
        end
        Fpd{paret_front_dim}(:,:,ip)=fbs;
        Xpd{paret_front_dim}(:,:,ip)=xbs;

        %% Plots - sub optimization results
        if logplot==1
        if nfun==0%4
            figure(hplotF)
                scatter3(fbs(:,1),fbs(:,2),fbs(:,3),50,fbs(:,4),'filled')
                colorbar
                grid on
                title(['F - ' MOptype ', pareto iter:' num2str(ipar)])
        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ft=f;
        %f=fis;
        if 1==0
            %feval(con,x,'p');
            %disp([0 1 f2; b f;1 0 f1])
            superf(f);
            plotsbestas=0;
            if plotsbestas
                figure
                grid on
                %scatter3(fis(:,1),fis(:,2),fis(:,3),40,fis(:,3),'filled')
                %ff=[f1; f2;f3];
                %bf=b*ff;
                np0=size(p0,1);
                hold on
                    scatter3(f(:,1),f(:,2),f(:,3),40,f(:,3),'filled')
                for i=1:np0
                    plot3(p0(i,1),p0(i,2),p0(i,3),'x')
                    text(p0(i,1),p0(i,2),p0(i,3),num2str(i))
                    %text(fi(i,1),fi(i,2),fi(i,3),num2str(i))
                end
                hold off
                figure
                scatter3(f(:,1),f(:,2),f(:,3),40,f(:,3),'filled')
            end
        end
    
    end
end