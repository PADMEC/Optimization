  function varargout=RBM_NLfesol(area,props,fext,glb,Zd,tipo)
%
% Realisa Análise via o MEF
% [u,sig,esf,vol,du,Lambd,Dvib,freq]

    ang = props(1,:);
	comp = props(2,:);
	els = props(3,:); 
	ro = props(4,:);
    
    volem=area.*comp;
    TotReg=sum(volem);
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
F=fext*0;
E=comp'*0+E0;

en=0;
fext0=fext;ip=0;
max_iter=100;
while nl_tol
    iter=iter+1;
    R=fext-F;
    %plot_sig_vet(R,u,sigi,1)
    if tipo(1)>0 |tipo(2)>0

        %K0u=F, K0.u0=F0=F-R, K0.u0+[dF/du]Du-F=0, 
        %F0+dF-F=0, dF=-R, dK*du=R, 
        %
        %[dF/du]Du=dK*du=R
        %
        %[dK/du]
        % Matriz de rigidez tangente:
        %
        %[dkg, iglnr, bars_out] = keglb(area,comp,E,ang,glb);
        %
        % Deslocamentos:
        %
        du=u*0;
        %rM=rank(full(dkg(iglnr,iglnr)));
         %if rcond(full(dkg(iglnr,iglnr)))<1e-16
         %    con_tol=0;
         %    nl_tol=0;
%             figure, plot_sig_vet(F,u,sig,tp,bars_out,iglnr)
         %end
        %dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
        if iter==1
            [dkg, iglnr, bars_out] = keglb(area,comp,E,ang,glb);
            dKn=Zd(iglnr,:)'*dkg(iglnr,iglnr)*Zd(iglnr,:);
            dKn1=inv(dKn);
        end
        Rn=Zd(iglnr,:)'*R(iglnr);
        alf = dKn1*Rn;
        du(iglnr) = Zd(iglnr,:)*alf;
        nnu=find(isnan(du));
        du(nnu)=0;
        u=u+du;
    end
    if tipo(1)>0

        %
        % Energia
        %
        en=en+fext'*u;
        
        %
        % Deformacoes:
        %
        rstep=0;rst=1;
        while rstep<1
            es = trdeform( comp, ang, glb, u);
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
             if rstepl<1
                rstepl
                fext=fext*rstepl;
                u=u*rstepl;
                es = trdeform( comp, ang, glb, u);
                [sig,E]=strain2stress(mat_curv,es);
                fprintf('Escoamento, tensao maxima: %f \n',max(abs(sig)));
                PlastRegion=0;
                rstepl=rstepl*rst;
             end
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
        [F] = forcs(area, ang, glb, sig);
    

        %E=els;
        x=area;
        floc=sig./(-x'.*E*pi/4./comp'.^2);
        varargout={u,sig,en,vol,floc};
    end

    %flamb=1;
    u_flamb=[];Dvib=[];
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
        [freq, Dvib]=vibsol(props,glb,dkg,area,ro);

    end
    
    fprintf('step: %d, iter: %d, maximas, deformacao: %f, tensao: %f , R:%f, du:%f \n'...
        ,ip,iter, max(abs(es)),max(abs(sig)), norm(R), norm(du)/norm(u));
    
    %%%%%%%%%%%%%% Bars Out %%%%%%%%%%%%%
    if ~isempty(bars_out)
        figure, plot_sig_vet(F,u,sig,1,bars_out,iglnr)
        nl_tol=0;
    end

    du_tol=norm(du)/norm(u);
    Rn_tol=norm(R)/norm(fext0);
    if (du_tol<1e-5)&(Rn_tol<1e-5)
        
        if norm(fext)==norm(fext0)
            nl_tol=0;
        end
        
        fext=fext+fext0*rstepl/10;
        if norm(fext)>norm(fext0)
            fext=fext0;
        end
        
        ip=ip+1;
        Fs(:,ip)=fext;
        Us(:,ip)=u;
        Sigs(:,ip)=sig;
        iter=0;
        %figure, 
        %plot_sig_vet(R,u,sig,1)
    end
    
    if (norm(fext0-fext))/norm(fext0)<1
        %nl_tol=1;
    else
        nl_tol=0;
    end
    
    if iter>max_iter
        nl_tol=0;
    end
        
    
%     ms=mean(abs(sig));
%     Dmean_S=abs(mean_sig-ms);
%     if Dmean_S<ms/1e3
%         %nl_tol=0;
%     end
%     mean_sig=ms;
    
    es0=es;
end

varargout={Fs,Us,Sigs,en,vol,floc,u_flamb,Lambd,Dvib,freq};