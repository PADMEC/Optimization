function varargout=NLfesol(area,props,fext,glb,tipan,ROM,Zd)
%
% Realisa Análise Nao Linear via o MEF
% [u,sig,esf,vol,du,Lamb,D,freq]

global offline_stage grad_calculation u_initial restart_parm ...
        link echo Truss3D NLGEOM

        %K0u=F, K0.u0=F0=F-R, K0.u0+[dF/du]Du-F=0, 
        %F0+dF-F=0, dF=-R, dK*du=R, 
        %
        %[dF/du]Du=dK*du=R
        % Matriz de rigidez tangente:
        
    cosn = props{1};
	comp = props{2};
	els = props{3};
	ro = props{4};
    
    volem=area.*comp;
    TotReg=sum(volem);
    Ulength=0;

if nargin<7
    ROM=0;
end
if nargin<5
    tipan=[1 0 0];
end

if sum(tipan)<0
    tipan=abs(tipan);
    sens=1;
else
    sens=0;
end

global es0  mat_curv
is_linearMat=0;
nl_tol=1;mean_sig=0;
ngl=max(max(glb));
K0=zeros(ngl,ngl);
PlastRegion=0;
fext0=fext;

if NLGEOM
    fexti=fext/4;
    fext=fexti;
end

[sigi0,E0]=strain2stress(mat_curv,1e-11,els);

%Strain Transormation (T): strain=T*u
%Normal force Transormation (T): F=T*norm
if Truss3D
    Tu_strain = Tstrain3D( comp, cosn, glb);
%     Mesf_F = T_esf_F3D( cosn, glb);
else
    Tu_strain = Tstrain  ( comp, cosn, glb);
end
Mesf_F = T_esf_F  ( cosn, glb);

if (ROM>=1)
    %Tu_strainN = Tu_strain*Zd;
    alf=zeros(size(Zd,2),1);
    K0n=zeros(size(Zd,2));
end

if NLGEOM, 
    flex=[];
    esf=comp*0;
    k0 = zeros(ngl,ngl);
    is_linear=0;
end

%grad_calculation=0;
if grad_calculation
    %restart_parm = [is_nonlinear, factor, rstepl];
    is_nonlinear=restart_parm(1);
    if is_nonlinear
        is_linear=0;
        factor=restart_parm(2);
        rstepl=restart_parm(3);
        fext=fext/factor;
        ip=0;
        u=u_initial;
        es = Tu_strain*u;
        [sig,E]=strain2stress(mat_curv,es,els');
        esf = sig.*area';
        F=Mesf_F*esf;
        en=fext'*u;
        R=u*0;
        if echo
            fprintf('- starting step: %d, maximas, fext: %f, F: %f, tensao: %f \n'...
             ,ip, max(abs(fext)), max(abs(F)),max(abs(sig)));      
        end
    end
else is_nonlinear=0;
end
if ~is_nonlinear
    ip=0;
    u=zeros(ngl,1);
    F=fext*0;
    E=comp'*0+E0';
%     R=u;
%     dE_elm=E;
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
%     K_atual=0;
    K_atual=1;
else
    %K_atual=0;
    K_atual=1;
end

while nl_tol
%% Start Analysis
    iter=iter+1;
    R=fext-F;
    if iter==1
        %R=fext*1e-12;
    end
    %plot_sig_vet(R,u,sigi,1)
    if tipan(1)>0 ||tipan(2)>0

        % Deslocamentos:
        %
        du=u*0;
         %rM=rank(full(dkg(iglnr,iglnr)));
         %if rcond(full(dkg(iglnr,iglnr)))<1e-16
         %    con_tol=0;
         %    nl_tol=0;
%             figure, plot_sig_vet(F,u,sig,tp,bars_out,iglnr)
         %end
        if (ROM>0)
            %% ON_LINE STAGE
            if iter==1&&ip==0
                % 'u0' calculate
                
                %ui=~Zd.ai
                %u0+du=ui
                %u0+Zd.dai=ui
                %R=fext-K.u0
                %Kt.du=R 
                %Kt.Zd.dai=R
                
                [dkg, iglnr, bars_out] = keglb(area,comp,dE_elm,cosn,glb,Kelem_new);
                if NLGEOM,
                    K0=dkg+k0;
                else
                    K0=dkg;
                end
                dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
                %dRu=dkg(iglnr,iglnr)*(u(iglnr));
                Rn=Zd(iglnr,:)'*(R(iglnr));
                dalf = dKn\Rn;
                du(iglnr) = Zd(iglnr,:)*dalf;
                
                if ~K_atual
                    %dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
                    dKn1=inv(dKn);
                end
                K0n=dKn;
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
                    
%                         K0(gbs,gbs)=K0(gbs,gbs)+dkg(gbs,gbs);
%                         dkg=K0;
%                         dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);

%                         K0=K0+dKn;
%                         dKn=K0;

                        % Full K Calculation TEST
%                         [dkg1, iglnr1, bars_out1,gbs1] = keglb(area,comp,E,cosn,glb,1:2060);
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
                    
                    %(K+dK)*du=(Ku)+(Kdu)+dK*u+dK*du=~ R0 + Kdu + dR
                    % = Ri, Kdu = Ri - 
                    %dRu=dRu+prevR+dkg(iglnr,iglnr)*u(iglnr);
                    
                    %F0+Kt*du=Fe, Kt*du=R
                    %ROM -> du=Z.da
                    %Kt.da=R
                    
                        
                    [dkg, iglnr, ~,gbs] = keglb(area,comp,dE_elm,cosn,glb,Kelem_new);
                    if NLGEOM,
                        dkg=dkg+k0;
                    end
                    dpKn=Zd(gbs,:)'*dkg(gbs,gbs)*Zd(gbs,:);
                    dKn = K0n + dpKn;
                    
                    %dRu=dkg(iglnr,iglnr)*u(iglnr);
                    Rn=Zd(iglnr,:)'*R(iglnr);
                    
                    dalf = dKn\Rn;
                    
                    K0n = dKn;
                    
                    
                else
                    if iter==1
                        [dkg, iglnr, ~,gbs] = keglb(area,comp,dE_elm,cosn,glb,Kelem_new);
                        if NLGEOM,
                            dkg=dkg+k0;
                        end
                        K0(gbs,gbs)=K0(gbs,gbs)+dkg(gbs,gbs);
                        dkg=K0;
                        dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
                        dKn1=inv(dKn);
                    end
                    %dRu=dkg(iglnr,iglnr)*u(iglnr);
                    Rn=Zd(iglnr,:)'*(R(iglnr));
                    dalf = dKn1*Rn;
                end
                %du(iglnr) = Zdu(iglnr,:)*alf;
                du(iglnr) = Zd(iglnr,:)*dalf;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ROM Test
%             du2=du*0;
%             du2(iglnr) = dkg(iglnr,iglnr)\R(iglnr);
%             disp(norm(du-du2)/norm(du2));
%             alf0=alf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        elseif (ROM<=0)
            %% Single FEM or OFF_LINE STAGE
            [dkg, iglnr, ~,gbs] = keglb(area,comp,dE_elm,cosn,glb,Kelem_new);
            if NLGEOM,
                dkg=dkg+k0;
            end
            K0(gbs,gbs)=K0(gbs,gbs)+dkg(gbs,gbs);
            dkg=K0;
            du(iglnr) = dkg(iglnr,iglnr)\R(iglnr);
            
            % Full K Calculation TEST
            %[dkg1, iglnr1, bars_out1,gbs1] = keglb(area,comp,E,cosn,glb,1:2060);
            %norm(dkg-dkg1)

            if offline_stage
                if iter==1&&ip==0
                    K0=dkg;
                end
            end
            
        end
        nnu=(isnan(du));
        du(nnu)=0;
        u=u+du;
        
        if (ROM>0)
            alf = alf + dalf;
        end
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% pos-processamento %%%%
%             es = Tu_strain*u;
%             [sig,E]=strain2stress(mat_curv,es,els);
%             esf = sig.*area';
%             F=Mesf_F*esf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if tipan(1)>0

        %
        % Energia
        %
        en=en+R'*du;
        
        %
        % Deformacoes:
        %
        rstep=1;rst=1;
        
        %if (ROM==1)
        %    es = Tu_strainN*alf;
        %else
            es = Tu_strain*u;
        %end
        while rstep<1
            %es = trdeform( comp, cosn, glb, u);
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
        
        %
        % Tensões:
        %
        prevE=E;
        [sig,E]=strain2stress(mat_curv,es,els');
        
        if (iter==1)&&(ip==0)
            lin_sig=es.*E0';
            rstepl=mat_curv(1,2)*1.0001/max(abs(lin_sig));
             if rstepl<1
                 % PLASTIFICATION, BACK TO THE YELD POINT
                is_linearMat=0;
                %rstepl
                fext=fext*rstepl;
                u=u*rstepl;
                %es = trdeform( comp, cosn, glb, u);
                es = Tu_strain*u;
                [sig,E]=strain2stress(mat_curv,es,els');
                %fprintf('Escoamento, tensao maxima: %f , k: %f\n',max(abs(sig)),rstepl);
                %PlastRegion=0;
                rstepl=rstepl*rst;
             else
                 is_linearMat=1;
             end
        end
        
        %if is_linearMat~=1
            % Atualizacao de K onde a tensao mudou na curva
            dE_elm=E-prevE;
            Kelem_new=find(dE_elm);
        %end
        
        %
        % Esforços:
        %
        esf = sig.*area';%sigma(area,comp,esf);
        
        %
        %K(u)*u=
        %elastic: k0.u=f
        %plastic: F(u)=k0u0+dk.du, (dk=dF/du)
        %[F] = forcs(area, cosn, glb, sig);
        F=Mesf_F*esf;
        
        
        if NLGEOM
            k0 = kgeom(esf,comp,cosn,glb);
            F = F + k0*u;
        end
    

        %E=els;
        x=area;
        floc=sig./(-x'.*E*pi/4./comp'.^2);
        
    end

    %flamb=1;
    V=[];D=[];
    Lamb=zeros(length(fext),1);
    Lamb=zeros(min([length(fext),6]),1);
    freq=zeros(length(fext),1);

    if tipan(2)>0
    %floc=area.*els*pi/4./comp.^2;
    %if nargout>7
    %    [Lamb, du,uf,sigf,esff]=flambsol(esf,props,glb,kg,area,fext,u);
    %    varargout={u,sig,en,vol,floc,du,Lamb,uf,sigf,esff};
        [Lamb, V]=flambsol(esf,props,glb,dkg);
    end
    
    if tipan(3)>0
        [freq, D]=vibsol(props,glb,dkg,area,ro);
    end
    
    % COMPUTE TOLERANCES
    du_tol=norm(du)/norm(u);
    %Rn_tol=norm(R)/norm(fext0);
    
    
    if (ROM>0)
        Rn=Zd'*((F-fext));
        Rn_tol = max(abs(Rn))/max(abs(Zd'*fext));
    else
        Rn_tol = max(abs(F-fext))/max(abs(fext));
    end
    
    %Limite U Length
    if Ulength==0
        Ulength=norm(u);
        du=0;
    end
    
    if NLGEOM
        flex=[flex; [norm(u) norm(F) norm(fext)]];
    end
    
    
    %% Convergence test
    is_linear=is_linearMat&&~NLGEOM;
    if (du_tol<1e-6)&&(Rn_tol<1e-5)||(iter>max_iter)||is_linear
        if echo
        fprintf('- step: %d, iter: %d, maximas, F: %f, tensao: %f , R:%g, du:%g \n'...
          ,ip,iter, max(abs(fext)),max(abs(sig)), Rn_tol, du_tol);      
        end 

        %Step CONVERGED -> next step
        ip=ip+1;
        Fs(:,ip)=fext;
        Us(:,ip)=u;
        Sigs(:,ip)=sig;
        Lambs(:,ip)=Lamb;
        iter=0;
        T(ip)=cputime;
        
            Ulength=norm(u);
            du=0;
            
        if offline_stage
            ZK{ip}=sparse(dkg - K0);
            nnzKs(ip)=nnz(ZK{ip});
        end
        %disp('-')
        
        fext_old=fext;
        if norm(fext)>=norm(fext0)
            nl_tol=0;
        else

            % Increment computation
            if is_linearMat&&NLGEOM
                fext=fext+fexti;
            elseif K_atual==1
                fext=fext_old+fext0*rstepl/10;
            else
                fext=fext_old+fext0*rstepl/15;
            end
            if norm(fext)>norm(fext0)
                fext=fext0;
            end

            if max(abs(sig))>=mat_curv(end,2)*.999
                % ROMPEU
                nl_tol=0;
                disp(['Failure, sig = ' num2str(max(abs(sig)))])
            end
            
        end
            
        %figure, 
        %plot_sig_vet(area,F,u,sig,tp,bars_out,iglnr)
        
%     elseif Rn_tol>2
        %fext=fext_old+fext0*rstepl/10/2;
        %u0
        %sig0
        %F0
    else 
        stepsize=norm(du)/Ulength;
        if stepsize>1
            fprintf('Diverging at step %d, iter %d, decreassing step size %f. \n',...
                ip, iter, stepsize)
            if ip==0
                fext=(fext)/2;
                u=u*0;
                F=F*0;
                sig=sig*0;
            elseif ip==9999
                fext=(Fs(:,ip)+fext)/2;
                u=Us(:,ip);
                F=Fs(:,ip);
                sig=Sigs(:,ip);
                Lamb=Lambs(:,ip);
%                 iter=0;
            end
%             iter=0;
            
            if NLGEOM
                esf = sig.*area';%sigma(area,comp,esf);
                k0 = kgeom(esf,comp,cosn,glb);
            end
        end
    end
    
    if iter>max_iter
        nl_tol=0;
    end
    
    es0=es;
end
u_initial=u;
fext=fext_old;
%varargout={u,sig,en,vol,floc};
if echo
fprintf(' -- step: %d, iter: %d, max F: %g, tensao: %g , R:%g, du:%g \n'...
        ,ip,iter, max(abs(fext)),max(abs(sig)), Rn_tol, du_tol);
end   

%% Post-Processor

%
% Volume:
vol = sum(area.*comp.*ro);%volume(area,comp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Restarting parameters
mFext=max(abs(fext0));
mFs=max(abs(fext));
factor=(mFext/mFs);
is_nonlinear=~is_linear; %?????????????????
% restart_parm = [is_nonlinear, factor, rstepl];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify the end of the nonlinear analysis
if mFs<mFext
    % Case failure before the required load
    % majour the outputs to create gradients
    ip=ip+1;
    Fs(:,ip)=fext0;
    Us(:,ip)=u*factor;
    Sigs(:,ip)=sig*factor;
    T(ip)=cputime;
    en=en*factor;
    fprintf('step: %d, failure, maximas, F: %g, tensao: %g , R:%g, du:%g \n'...
            ,ip, max(abs(fext)),max(abs(sig)), Rn_tol, du_tol);
end
% esf = Sigs(:,ip).*area';
if ROM<0
    %OFF_LINE STAGE
    varargout={Fs,Us,Sigs,en,vol,T,floc,ZK};
else
    varargout={Fs,Us,Sigs,en,vol,T,floc,V,-Lamb,D,freq};
end

global RandVarProp rv_data_comp var_type

if sens
%% Sensibility analysis
%   D(fi+fe)/Dx = 0 = dFi0/du * du/dx + dFe/dx 
%   dFi/dx + K * du/dx = - dF/dx
% x = rand var: 
% E -> du/dE = - K^-1*(dFi/dE) = -K^-1*(dK*u)
    
    if var_type==1
    % Design Variables sensibility - Cross-section Areas
        ndv=max(link(:,1));
        Var_kinds = ones(1,ndv);
        [dFi,dFe,dvol]=grad_FdA(link,comp,ro,ngl,sig,Mesf_F,1:ndv);
        dKt(1,1,1:ndv) = 1:ndv;

    else
        % Random Variables sensibility nd d2f/dudx
        nrv=length(rv_data_comp);
        dvol=zeros(nrv,1);
        % dFi/dx init
        dFi=zeros(ngl,nrv);
        % dFe/dx init
        dFe=zeros(ngl,nrv);
        
        i_var_mat = 0; % n of  material variables
        dKt = zeros(ngl,ngl,nrv); %material gradient (only for material variables)
        
        Var_kinds = RandVarProp(1,:);
        
        vargs={area,props,fext,glb,tipan};
        for iv=1:nrv
            %dxdx = zeros(1,nrv);
            %dxdx(iv)=1;
            components=rv_data_comp{iv};
            RVwich = RandVarProp(2,iv);
            RVkind = Var_kinds(iv);
            
            switch RVkind
            case 1
                % Design Random Variables
                if RVwich==5||RVwich==51
                    % Design Variable Factor!!!
                    %F(a) = F(a*fat)
                    %dF/dfat = dF/da*da/dfat = dF/da*a
                    i_bar=find(link(:,1)==components,1);
                    fat=area(i_bar);
                else
                    fat=1;
                end
                [dF,dFei,dvol(iv)]=grad_FdA(link,comp,ro,ngl,sig,Mesf_F,components);
                dFi(:,iv)=dF*fat;
                dFe(:,iv)=dFei*fat;
                dvol(iv)=dvol(iv)*fat;
                dKt(1,1,iv)=components;
                
            case 2
                % Material Random Variables
                iwich = 2;%element component
                if RVwich==1||RVwich==51
                    % Related to Design Variable!!!
                    components=find(link(:,1)==components);
                end
                
                if RVwich==5||RVwich==51
                    % Design Variable Factor!!!
                    %F(a) = F(a*fat)
                    %dF/dfat = dF/da*da/dfat = dF/da*a
                    fat=E(components(1));
                else
                    fat=1;
                end
                
                vargs{RVkind}{3}=E*0;
                updated_varg = RV_update(1,RVkind,iwich,...
                    components,vargs{RVkind});
                dE_elm = updated_varg{3};
                i_var_mat = i_var_mat + 1;
                
                [dKt(:,:,iv)] = keglb(area,comp,dE_elm,cosn,glb,0);
                if NLGEOM
                    %Esf = es*E = [d].u.E
                    %    es = Tu_strain*u;
                    %dEsf/dE =  [Tu_strain].u.dE_elm
                    %dEsf/du =  [Tu_strain].du.E_elm
                    dEsf =  (Tu_strain*u).*dE_elm;
                    k0 = kgeom(dEsf,comp,cosn,glb);
                    dKt(:,:,iv)=dKt(:,:,iv)+k0;
                end
                
                dFi(:,iv) = dKt(:,:,iv)*u*fat;
                
                %display('... Material Random Variable - Df/dx implementation needed (NLfesol)')

            case 3
                % Load Random Variables
                %   K * du/dfi =  dF/dfi
                % dF/df
                sgnF=sign(components);
                dFe(abs(components),iv)=sgnF;

            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%  GENERAL SENSITIVITY EQUATION  %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %d(Fi=Fe)/dx, (dFi/dx) + (dFi/du)*(du/dx) = dFe/dx
    % K * du/dx = dFe/dx - dFi/dx
    dF = dFe - dFi;
    % du/dx computing
    du = K0\dF;

    % d(K * du/dx = dFe/dx - dFi/dx)/dy
    % dK/dy*du/dx + K*ddu/dxdy = ddFe/dxdy - ddFi/dxdy
    
    [dsig,dF,den,desf]=grad_postpoc(du, u, sig, fext, E, area, link, Tu_strain, Mesf_F, Var_kinds);
    % Critic Load Gradient
    if tipan(2)>0
        [dLamb] = grad_flamb(Lamb,V,desf,K0,props,E,glb,link,Var_kinds,dKt);
    else 

        dLamb=zeros(length(den),min([length(fext),6]));
    end
    varargout={varargout{:},dF,du,dsig,den,dvol,-dLamb};
    
    %figure, plot_sig_vet(area,F,du(:,1),dsig(:,1),1)
    %figure, plot_sig_vet(area,F,u,sig,1)
    
    if (ROM>=1)
    elseif (ROM<=0)
        %dkg=K0;
    end

end

end
    