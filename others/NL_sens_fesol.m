  function varargout=NLfesol(area,props,fext,glb,tipo,ROM,Zd,Zkn,Cs,xi)
%
% Realisa Análise Nao Linear via o MEF
% [u,sig,esf,vol,du,Lambd,D,freq]

global offline_stage grad_calculation u_initial restart_parm

        %K0u=F, K0.u0=F0=F-R, K0.u0+[dF/du]Du-F=0, 
        %F0+dF-F=0, dF=-R, dK*du=R, 
        %
        %[dF/du]Du=dK*du=R
        % Matriz de rigidez tangente:
        
    ang = props(1,:);
	comp = props(2,:);
	els = props(3,:); 
	ro = props(4,:);
    
    volem=area.*comp;
    TotReg=sum(volem);

    
    
if nargin<7
    ROM=0;
end
if nargin<5
    tipo=[1 0 0];
end

global es0  mat_curv

nl_tol=1;mean_sig=0;
ngl=max(max(glb));
K0=zeros(ngl,ngl);
PlastRegion=0;
fext0=fext;

[sigi0,E0]=strain2stress(mat_curv,1e-11,els);

%Strain Transormation (T): strain=T*u
Tu_strain = Tstrain( comp, ang, glb);
%Normal force Transormation (T): F=T*norm
Mesf_F = T_esf_F( ang, glb);
%grad_calculation=0;
if grad_calculation
    %restart_parm = [is_nonlinear, factor, rstepl];
    is_nonlinear=restart_parm(1);
    if is_nonlinear
        is_linear=0;
        factor=restart_parm(2);
        rstepl=restart_parm(3);
        fext=fext/factor;
        ip=1;
        u=u_initial;
        es = Tu_strain*u;
        [sig,E]=strain2stress(mat_curv,es,els');
        esf = sig.*area';
        F=Mesf_F*esf;
        en=fext'*u;
        R=u*0;
         fprintf('- starting step: %d, maximas, fext: %f, F: %f, tensao: %f \n'...
             ,ip, max(abs(fext)), max(abs(F)),max(abs(sig)));      
    end
else is_nonlinear=0;
end
if ~is_nonlinear
    ip=0;
    u=zeros(ngl,1);
    R=u;
    F=fext*0;
    E=comp'*0+E0';
    dE_elm=E;
    en=0;
end

iter=0;
Kelem_new=1:length(E);
dE_elm=E;
max_iter=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose NL solver
if ROM>0
    max_iter=100;
    %Atualizar K - Matriz de Rigidez tg
    K_atual=0;
    K_atual=1;
else
    %K_atual=0;
    K_atual=1;
end

%% Start Analysis
while nl_tol
    iter=iter+1;
    prevR=R;
    R=fext-F;
    if iter==1
        %R=fext*1e-12;
    end
    %plot_sig_vet(R,u,sigi,1)
    if tipo(1)>0 ||tipo(2)>0

        % Deslocamentos:
        %
        du=u*0;
        %rM=rank(full(dkg(iglnr,iglnr)));
         %if rcond(full(dkg(iglnr,iglnr)))<1e-16
         %    con_tol=0;
         %    nl_tol=0;
%             figure, plot_sig_vet(F,u,sig,tp,bars_out,iglnr)
         %end
        if (ROM==1)
            %% ON_LINE STAGE
            if iter==1&&ip==0
                % 'u0' calculate
                
                %u=~Zd.a!
                %u0+du=u!
                %R=fext-Ku0
                %Kt.du=R! ->Kt.(u-u0)=R
                %Kt.u-Kt.u0=R
                %Kt.Zd.a=~R+Kt.u0
                %Rt=R+Kt.u0! (Ku=Fi, R=F0-Ku0)
                %Kt.Zd.a=Rt
                %Zd'.Kt.Zd.a=Zd'.Rt
                %u0+du=Zd.a
                
                [dkg, iglnr, bars_out] = keglb(area,comp,dE_elm,ang,glb,Kelem_new);
                K0=dkg;
                dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
                dRu=dkg(iglnr,iglnr)*u(iglnr);
                Rn=Zd(iglnr,:)'*(R(iglnr)+dRu);
                alf = dKn\Rn;
                %alf = alf*0;alf(3) =1;
                du(iglnr) = Zd(iglnr,:)*alf - u(iglnr);
                
                if ~K_atual
                    %dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
                    dKn1=inv(dKn);
                end
                %K0=dKn;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            
%                 dKn*(alf0+alf) = Rn;
%                 alf = dKn1*(Rn+dKn*alf0);
%                 unew(iglnr) = Zd(iglnr,:)*alf;
%                 du = unew-u;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                % 'du' calculate
                
                %Rn=Zdu(iglnr,:)'*R(iglnr);
                
                if K_atual
                        
                        [dkg, iglnr, bars_out,gbs] = keglb(area,comp,dE_elm,ang,glb,Kelem_new);
                        K0(gbs,gbs)=K0(gbs,gbs)+dkg(gbs,gbs);
                        dkg=K0;
                        dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
%                         K0=K0+dKn;
%                         dKn=K0;

                        % Full K Calculation TEST
%                         [dkg1, iglnr1, bars_out1,gbs1] = keglb(area,comp,E,ang,glb,1:2060);
%                         dKn1=Zd(iglnr,:)'*dkg1(iglnr,iglnr)*Zd(iglnr,:);
%                         norm(full(dKn-dKn1))
                    
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         % K aproximation (g weigth)
%                         Xi=[max(abs(F)) xi];
%                         b=gWeigth(Cs',Xi);
%                         nk=length(Zkn);ddKn=Zkn{1}*0;
%                         for ik=1:nk
%                             ddKn=ddKn+Zkn{ik}*b(ik);
%                         end
%                         dKn=K0+ddKn;
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     dRu1=dkg1(iglnr,iglnr)*u(iglnr);
%                     Rn1=Zd(iglnr,:)'*(R(iglnr)+dRu1);
                    
                    %(K+dK)*(u+du)=(Ku)+(Kdu)+dK*u+dK*du=dRu+R+
                    %dRu=dRu+prevR+dkg(iglnr,iglnr)*u(iglnr);
                    dRu=dkg(iglnr,iglnr)*u(iglnr);
                    Rn=Zd(iglnr,:)'*(R(iglnr)+dRu);
                    
                    %norm(full(dRu-dRu1))
                    %norm(full(Rn-Rn1))
                    
                    alf = dKn\Rn;
                else
                    if iter==1
                        [dkg, iglnr, bars_out,gbs] = keglb(area,comp,dE_elm,ang,glb,Kelem_new);
                        K0(gbs,gbs)=K0(gbs,gbs)+dkg(gbs,gbs);
                        dkg=K0;
                        dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
                        dKn1=inv(dKn);
                    end
                    dRu=dkg(iglnr,iglnr)*u(iglnr);
                    Rn=Zd(iglnr,:)'*(R(iglnr)+dRu);
                    alf = dKn1*Rn;
                end
                %du(iglnr) = Zdu(iglnr,:)*alf;
                du(iglnr) = Zd(iglnr,:)*alf - u(iglnr);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             du2=du*0;
%             du2(iglnr) = dkg(iglnr,iglnr)\R(iglnr);
%             disp(norm(du-du2)/norm(du2));
%             alf0=alf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        elseif (ROM<=0)
            %% Single FEM or OFF_LINE STAGE
            [dkg, iglnr, bars_out,gbs] = keglb(area,comp,dE_elm,ang,glb,Kelem_new);
            K0(gbs,gbs)=K0(gbs,gbs)+dkg(gbs,gbs);
            dkg=K0;
            du(iglnr) = dkg(iglnr,iglnr)\R(iglnr);
            
            % Full K Calculation TEST
            %[dkg1, iglnr1, bars_out1,gbs1] = keglb(area,comp,E,ang,glb,1:2060);
            %norm(dkg-dkg1)
            
            %TEST u SPACE
                % Fint=Fext
                % F0+du*F'=Fext
                % dK(ui-u0)=R
                % dK*ui=R+dK*u0
                % dK*ui=R+dR
                
%                 dRu=dkg(iglnr,iglnr)*u(iglnr);
%                 Rn=(R(iglnr)+dRu);
%                 ui = dkg(iglnr,iglnr)\Rn;
%                 du(iglnr) = ui - u(iglnr);

            if offline_stage
                if iter==1&&ip==0
                    K0=dkg;
                end
            end
                
            
%             %[x,flag,Rr,it,Rn]=method(A,b,tol,maxit,M1,M2,x0)
%             A=dkg(iglnr,iglnr);
%             b=R(iglnr); x0=du(iglnr);
%             [x] = pcg(A,b,x0); %pcg
%             %[L,U] = lu(A); x = Un(Lnb); %For sparse A
%             du(iglnr)=x;
        end
        nnu=(isnan(du));
        du(nnu)=0;
        u=u+du;
        
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% pos-processamento %%%%
%             es = Tu_strain*u;
%             [sig,E]=strain2stress(mat_curv,es,els);
%             esf = sig.*area';
%             F=Mesf_F*esf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if tipo(1)>0

        %
        % Energia
        %
        en=en+fext'*du;
        
        %
        % Deformacoes:
        %
        rstep=0;rst=1;
        while rstep<1
            %es = trdeform( comp, ang, glb, u);
            es = Tu_strain*u;
            iPLastEel=(abs(es)>mat_curv(1));
            PlastRegion1=sum(volem(iPLastEel));
%            if PlastRegion1==0
                rstep=2;
%            else
%                 rstep=(PlastRegion+.4*TotReg)/PlastRegion1;
%                 if rstep<1
%                     rstep
%                     fext=fext*rstep;
%                     u=u*rstep;
%                     rst=rst*rstep;
%                 end
%            end
        end
        PlastRegion=PlastRegion1;
        %
        % Tensões:
        %
        prevE=E;
        [sig,E]=strain2stress(mat_curv,es,els');
        if (iter==1)&&(ip==0)
            lin_sig=es.*E0';
            rstepl=mat_curv(1,2)/max(abs(lin_sig))*.999;
             if rstepl<1
                 % PLASTIFICATION, BACK TO THE YELD POINT
                is_linear=0;
                %rstepl
                fext=fext*rstepl;
                u=u*rstepl;
                %es = trdeform( comp, ang, glb, u);
                es = Tu_strain*u;
                [sig,E]=strain2stress(mat_curv,es,els');
                %fprintf('Escoamento, tensao maxima: %f , k: %f\n',max(abs(sig)),rstepl);
                PlastRegion=0;
                rstepl=rstepl*rst;
             else
                 is_linear=1;
             end
        end
        
        % Atualizacao de K onde a tensao mudou na curva
        dE_elm=E-prevE;
        Kelem_new=find(dE_elm);
        
        %
        % Esforços:
        %
        esf = sig.*area';%sigma(area,comp,esf);
        %
        % Volume:
        %
        vol = sum(area.*comp.*ro);%volume(area,comp);
        
        %
        %K(u)*u=
        %elastic: k0.u=f
        %plastic: F(u)=k0u0+dk.du, (dk=dF/du)
        %[F] = forcs(area, ang, glb, sig);
        F=Mesf_F*esf;
    

        %E=els;
        x=area;
        floc=sig./(-x'.*E*pi/4./comp'.^2);
        
    end

    %flamb=1;
    u_flamb=[];D=[];
    Lambd=zeros(1,length(fext));
    freq=zeros(1,length(fext));

    if tipo(2)>0
    %floc=area.*els*pi/4./comp.^2;
    %if nargout>7
    %    [Lambd, du,uf,sigf,esff]=flambsol(esf,props,glb,kg,area,fext,u);
    %    varargout={u,sig,en,vol,floc,du,Lambd,uf,sigf,esff};
        [Lambd, u_flamb]=flambsol(esf,props,glb,dkg);
    end
    
    if tipo(3)>0
        [freq, D]=vibsol(props,glb,dkg,area,ro);
    end
    
    du_tol=norm(du)/norm(u);
    %Rn_tol=norm(R)/norm(fext0);
    Rn_tol=max(abs(F-fext0))/max(abs(fext0));
    
    %% Convergence test
    if (du_tol<1e-6)||(Rn_tol<1e-4)||(iter==max_iter)||is_linear
%       fprintf('- step: %d, iter: %d, maximas, F: %f, tensao: %f , R:%g, du:%g \n'...
%           ,ip,iter, max(abs(fext)),max(abs(sig)), Rn_tol, du_tol);      
        
        %%%%%%%%%%%%%% Bars Out %%%%%%%%%%%%%
%         if ~isempty(bars_out)
%             figure, plot_sig_vet(F,u,sig,1,bars_out,iglnr)
%             nl_tol=0;
%         end

        %Step CONVERGED -> next step
        ip=ip+1;
        Fs(:,ip)=fext;
        Us(:,ip)=u;
        Sigs(:,ip)=sig;
        iter=0;
        T(ip)=cputime;

        if offline_stage
            ZK{ip}=sparse(dkg - K0);
            nnzKs(ip)=nnz(ZK{ip});
        end
        %disp('-')
        
        if norm(fext)>=norm(fext0)
            nl_tol=0;
        end
        fext_old=fext;
        
        if K_atual==1
            fext=fext_old+fext0*rstepl/15;
        else
            fext=fext_old+fext0*rstepl/20;
        end
        if norm(fext)>norm(fext0)
            fext=fext0;
        end

        if max(abs(sig))>=mat_curv(end,2)*.999
            nl_tol=0;
            disp('Failure')
        end
            
        %figure, 
        %plot_sig_vet(R,u,sig,1)
    elseif Rn_tol>2
        %fext=fext_old+fext0*rstepl/10/2;
        %u0
        %sig0
        %F0
        
    end
    
%     if iter>max_iter
%         nl_tol=0;
%     end
        
    
%     ms=mean(abs(sig));
%     Dmean_S=abs(mean_sig-ms);
%     if Dmean_S<ms/1e3
%         %nl_tol=0;
%     end
%     mean_sig=ms;
    
    es0=es;
end
u_initial=u;
fext=fext_old;
%varargout={u,sig,en,vol,floc};
fprintf(' -- step: %d, iter: %d, maximas, F: %g, tensao: %g , R:%g, du:%g \n'...
        ,ip,iter, max(abs(fext)),max(abs(sig)), Rn_tol, du_tol);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set Restarting parameters
        mFext=max(abs(fext0));
        mFs=max(abs(fext));
        factor=(mFext/mFs);
        is_nonlinear=~is_linear;
        restart_parm = [is_nonlinear, factor, rstepl];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Verify the end of the nonlinear analysis
        if mFs<mFext
            % Case failure before the required load
            % majour the outputs to create gradients
            ip=ip+1;
            %Fs(:,ip)=fext0;
            %Us(:,ip)=u*factor;
            %Sigs(:,ip)=sig*factor;
            %T(ip)=cputime;
            %en=en*factor;
            %fprintf('step: %d, failure, maximas, F: %g, tensao: %g , R:%g, du:%g \n'...
            %        ,ip, max(abs(fext)),max(abs(sig)), Rn_tol, du_tol);
        end

if ROM<0
    %OFF_LINE STAGE
    varargout={Fs,Us,Sigs,en,vol,T,floc,ZK};
else
    varargout={Fs,Us,Sigs,en,vol,T,floc,u_flamb,Lambd,D,freq};
end