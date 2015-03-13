function [dFi,dFe,dvol]=grad_FdA(link,comp,ro,ndof,sig,Mesf_F,components)
% Area cross section derivative

    % Design Variables sensibility - Cross-section Areas
    %
    %strain = Tu_strain*u;
    %[sig,E]=strain2stress(mat_curv,strain,els');
    % Fi(A) = M* esf = M*[ A ]*E*strain
    %dFi(A) = M*desf = M*[dA ]*E*strain
    %or:
    % Fi(A) = K*u, %dFi(A) = dK*u
    % [dKt] = keglb(dA,comp,E,cosn,glb,Kelem_new);
    
    ndv=length(components);
    nel=length(comp);

    % dFi/dx
    dFi=zeros(ndof,ndv);
    % dFe/dx
    dFe=zeros(ndof,ndv);
    dvol=zeros(ndv,1);
    dA=zeros(nel,ndv);
    rc=ro.*comp;
    for i=1:ndv
        ibv=link(:,1)==components(i);
        
        dA(:,i)=ibv;
        dvol(i,1) = sum(rc(ibv));
        
        % dFdx Computing
        desf = sig*0;
        desf(ibv) = sig(ibv);
        dFi(:,i)=Mesf_F*desf;
        %dKt = Mesf_F.*da'*E.*strain
    end
        