  function varargout=Lfesol(area,props,fext,glb,tipo,ROM,Zd,Zdu)
%
% Realisa Análise Nao Linear via o MEF
% [u,sig,esf,vol,du,Lambd,D,freq]

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
u=zeros(ngl,1);
iter=0;
PlastRegion=0;

[sigi0,E0]=strain2stress(mat_curv,1e-11);
Tu_strain = Tstrain( comp, ang, glb);
Mesf_F = T_esf_F( ang, glb);

F=fext*0;
E=comp'*0+E0;
en=0;
fext0=fext;ip=0;
fext_old=fext*0;
max_iter=10;
if ROM>0
    max_iter=100;
%     alf0=Zd(1,:)*0;
%     unew=zeros(ngl,1);
    K_atual=1;
else
    K_atual=1;
end
%while nl_tol
    iter=iter+1;
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
            %ON_LINE STAGE
            if iter==1&&ip==0
                % 'u0' calculate
                
                %u=~Zd.a
                %u0+du=u
                %Kt.du=R
                %Kt.u0+R=~Kt.Zd.a
                %Rt=Kt.u0+R (Ku=Fi, R=F0-Ku0)
                %Kt.Zd.a=Rt
                %Zd'.Kt.Zd.a=Zd'.Rt
                %u0+du=Zd.a
                
                [dkg, iglnr, bars_out] = keglb(area,comp,E,ang,glb);
                dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
                Rn=Zd(iglnr,:)'*(R(iglnr)+dkg(iglnr,iglnr)*u(iglnr));
                alf = dKn\Rn;
                du(iglnr) = Zd(iglnr,:)*alf - u(iglnr);
                
                if ~K_atual
                    dKn=Zdu(iglnr,:)'*dkg(iglnr,iglnr)*Zdu(iglnr,:);
                    dKn1=inv(dKn);
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                 alf = dKn1*(Rn+dKn*alf0);
%                 unew(iglnr) = Zd(iglnr,:)*alf;
%                 du = unew-u;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                % 'du' calculate
                
                Rn=Zdu(iglnr,:)'*R(iglnr);
                if K_atual
                    [dkg, iglnr, bars_out] = keglb(area,comp,E,ang,glb);
                    dKn=Zdu(iglnr,:)'*dkg(iglnr,iglnr)*Zdu(iglnr,:);
                    alf = dKn\Rn;
                else
                    if iter==1
                        [dkg, iglnr, bars_out] = keglb(area,comp,E,ang,glb);
                        dKn=Zdu(iglnr,:)'*dkg(iglnr,iglnr)*Zdu(iglnr,:);
                        dKn1=inv(dKn);
                    end
                    alf = dKn1*Rn;
                end
                du(iglnr) = Zdu(iglnr,:)*alf;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             du2=du*0;
%             du2(iglnr) = dkg(iglnr,iglnr)\R(iglnr);
%             disp(norm(du-du2)/norm(du2));
%             alf0=alf;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        elseif (ROM<=0)
            %single FEM or OFF_LINE STAGE
            [dkg, iglnr, bars_out] = keglb(area,comp,E,ang,glb);
            du(iglnr) = dkg(iglnr,iglnr)\R(iglnr);
            
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
        
        % 'Zdu' - OFF_LINE STAGE
        if (ROM<0)
            if ~((iter>1))
                Zdu=[Zdu du];
            end
        end
    end
    if tipo(1)>0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Criar Matrizes para pos-processamento %%%%
        %sig=Msig*u
        %F=MF*sig
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
            iPLastEel=find(abs(es)>mat_curv(1));
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
        [sig,E]=strain2stress(mat_curv,es);
        if (iter==1)&(ip==0)
            lin_sig=es*E0;
            rstepl=mat_curv(1,2)/max(abs(lin_sig));
             %if rstepl<1
                %rstepl
                fext=fext*rstepl;
                u=u*rstepl;
                %es = trdeform( comp, ang, glb, u);
                es = Tu_strain*u;
                [sig,E]=strain2stress(mat_curv,es);
                fprintf('Escoamento, tensao maxima: %f , k: %f\n',max(abs(sig)),rstepl);
                PlastRegion=0;
                rstepl=rstepl*rst;
             %end
        end
            
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
    Rn_tol=norm(R)/norm(fext0);
    fprintf('step: %d, iter: %d, maximas, F: %f, tensao: %f , R:%f, du:%f \n'...
        ,ip,iter, norm(fext),max(abs(sig)), Rn_tol, du_tol);
    
    %%%%%%%%%%%%%% Bars Out %%%%%%%%%%%%%
    if ~isempty(bars_out)
        %figure, plot_sig_vet(F,u,sig,1,bars_out,iglnr)
        %nl_tol=0;
    end

    if (du_tol<1e-6)||(Rn_tol<1e-4)||(iter==max_iter)
        
        if norm(fext)==norm(fext0)
            nl_tol=0;
        end
        fext_old=fext;
        
        if K_atual==1
            fext=fext_old+fext0*rstepl/3;
        else
            fext=fext_old+fext0*rstepl/18;
        end
        if norm(fext)>norm(fext0)
            fext=fext0;
        end
        if (du_tol<1e-3)||(Rn_tol<1e-2)
            %Step CONVERGED -> next step
            ip=ip+1;
            Fs(:,ip)=fext;
            Us(:,ip)=u;
            Sigs(:,ip)=sig;
        end
        iter=0;
        disp('-')
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
%end
%varargout={u,sig,en,vol,floc};
if ROM<0
    %OFF_LINE STAGE
    varargout={Fs,Us,Sigs,en,vol,Zdu};
else
    varargout={Fs,Us,Sigs,en,vol,floc,u_flamb,Lambd,D,freq};
end