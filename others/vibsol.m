function [freq,du]=vibsol(props,glb,kg,area,ro)

	comp = props{2};
    
    M = mmglb(area,comp,glb,ro);
    
    %kg*du+Lambd*mg*du==Indet
    %[D,W2]=eigs(kg,M);
    [D,W2]=eig(full(kg),full(M));
 
    W2=diag(W2);
    
    [Lambd,mdv]=sort(abs(real(W2)));
    %Lambd=full(Lambd);
    du=full((D(:,mdv)));
    
    freq=abs(Lambd).^.5;
