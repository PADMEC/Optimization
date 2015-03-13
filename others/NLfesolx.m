  function varargout=NLfesol(area,props,fext,glb,tipo)
%
% Realisa Análise via o MEF
% [u,sig,esf,vol,du,Lambd,D,freq]

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
figure,
while nl_tol
    iter=iter+1
        [sigi,E]=strain2stress(mat_curv,es0);
        [F] = forcs(area, ang, glb, sigi);
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
        [dkg, iglnr] = keglb(area,comp,E,ang,glb);
        
        %
        % Deslocamentos:
        %
        du=u*0;
        du(iglnr) = dkg(iglnr,iglnr)\R(iglnr);
        nnu=find(isnan(du));
        du(nnu)=0;
        u=u+du;
    end
    if tipo(1)>0

        %
        % Energia
        %
        en=fext'*u;
        
        %
        % Deformacoes:
        %
        rstep=0;
        while rstep<1
            es = trdeform( comp, ang, glb, u);
            iPLastEel=find(abs(es)>mat_curv(1));
            PlastRegion1=sum(volem(iPLastEel));
            if PlastRegion1==0
                rstep=2;
            else
                rstep=(PlastRegion+.1*TotReg)/PlastRegion1;
                if rstep<1
                    rstep
                    fext=fext*rstep;
                    u=u*rstep;
                end
            end
        end
        PlastRegion=PlastRegion1;
        %
        % Tensões:
        %
        [sig]=strain2stress(mat_curv,es);
        if iter==1
            lin_sig=es*E0;
            rstep=mat_curv(1,2)/max(abs(lin_sig));
             if rstep<1
                rstep
                %fext=fext*rstep;
                u=u*rstep;
                es = trdeform( comp, ang, glb, u);
                [sig]=strain2stress(mat_curv,es);
                fprintf('Escoamento, tenssao maxima: %f \n',max(abs(sig)));
                PlastRegion=0;
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

        %E=els;
        x=area;
        floc=sig./(-x'.*E*pi/4./comp'.^2);
        varargout={u,sig,en,vol,floc};
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
    if du_tol<1e-5
        nl_tol=0;
    end
    
    
    ms=mean(abs(sig));
    Dmean_S=abs(mean_sig-ms);
    if Dmean_S<ms/1e3
        %nl_tol=0;
    end
    mean_sig=ms;
    es0=es;
    fprintf('deformacao maxima: %f , R:%f, du:%f \n',max(abs(es)), norm(R), norm(du)/norm(u));
end

varargout={u,sig,en,vol,floc,u_flamb,Lambd,D,freq};