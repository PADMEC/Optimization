function [dLamb] = grad_flamb(Lamb,V,desf,Kg,props,E,glb,link,Var_kinds,dKtm)
% Comput the gradient of the critical load

    %%%%%%%%%%%%%%%%%%%%%%%%
    %Sensibility Stability
    %(Kt+Lambd*kg).V==0
    %(Kt+Lambd*kg).dV + (dKt+dLambd.kg+Lambd.dkg).V=0, [V'*(Kt+Lambd*kg)=0]
    %V'*(dKt+dLambd.kg+Lambd.dkg).V=0
    %dLambd.(V'*kg.V) + V'*(dKt+Lambd.dkg).V=0
    %dLambd * vKv = - V'*(dKt+Lambd.dkg).V=0
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    cosn  = props{1};
	comp  = props{2};
    ndv   = length(Var_kinds);
    nLamb = length(Lamb);
    nelem = length(comp);
    dLamb = zeros(ndv,nLamb);
    for il=1:nLamb
        % Computing: vKv
        vKv=(V(:,il)'*Kg*V(:,il));
        for ibv=1:ndv
        
            %%%%%%%%%%%%%%%%%%%%%%%%
            % dKt/dx
            if Var_kinds(ibv)==1
                dA= zeros(nelem,1);
                dA(link(:,1)==dKtm(1,1,ibv))=1;
                dKt = keglb(dA,comp,E,cosn,glb,0);
            elseif Var_kinds(ibv)==2
                %Material dKt computed
                dKt = dKtm(:,:,ibv);
            elseif Var_kinds(ibv)==3
                %External Force dKt = 0;
                dKt = Kg*0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%
            % dKg/dx
            dKg = kgeom(desf(:,ibv), comp, cosn, glb);

            %%%%%%%%%%%%%%%%%%%%%%%%
            % Computing: dLamb/dx
            dLamb(ibv,il)=V(:,il)'*(dKt+Lamb(il)*dKg)*V(:,il)/vKv*Lamb(il);
        end
    end
    