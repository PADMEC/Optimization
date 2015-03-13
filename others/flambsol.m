function [Lambd,du,uf,sigf,esff]=flambsol(esf,props,glb,kg,area,fext)
    
    cosn = props{1};
	comp = props{2};
	els = props{3}; 
    
    k0 = kgeom(esf,comp,cosn,glb);
    
    %kg*du+Lambd*kgeo*du==Indet
    %[V,D]=eigs(kg,-kgeo,1,'sm');
%     [V,D,erflag]=InvIter(kg,-k0);
%     if erflag
        [V,D]=eig(full(kg),full(-k0));
%     end
    
    
    D=diag(D);
    [Lambd,mdf]=sort(abs(real(D)));
    if (length(mdf)>5)
        Lambd(6:end)=[];
        mdf(6:end)=[];
    end
%     Lambd=Lambd(1);
%     mdf=mdf(1);
    %Lambd=full(Lambd);
    du=full((V(:,mdf)));
    if nargout>2
        if nargin>5
        dkg=kg+k0;
        u2=dkg\fext;
        uf=u2;
        %uf=u+u2;
        esff = tresf(area,comp,els,ang,glb,uf);               
        sigf = esff'./area';
        else
            warning('mais input: flambsol...')
        end
    end
    
    %Sensibility
    %(Kt+Lambd*kg).u==0 (q.q. u)
    %dKt.u+Kt.du+dLambd.kg.u+Lambd.dkg.u+Lambd.kg.du==Indet
    %(Kt + Lambd.kg).du + dLambd.kg.u = - (dKt.u + Lambd.dkg.u)
    