function [modosFLAMB,lambFLAMB,erflag] = InvIter(kg,kgg)

% Subrotina utilizada para calcular o menor autovalor e seu respectivo
% autovetor utilizando o metodo da iteracao inversa
% (kg + Lamb*kgg)*V=[0]
% Lamb = -v*inv(K)v/(vKgv)
lamb1 = 1;
v = ones(size(kg,1),1);
wb = kgg*v;

toler = 1e-8;
niter = 400;

invkg = inv(kg);

i = 1;
while i <= niter
    v = invkg*wb;
    dot1 = v'*wb;
    wb = kgg*v;
    dot2 = v'*wb;
    lamb2 = -dot1/dot2;
    
    Erro = abs(lamb2 - lamb1)/(1 + abs(lamb2));
    aux = sqrt(abs(dot2));
    wb = wb/aux;    
    lamb1 = lamb2 ;       
    if Erro <= toler
       break
   end
   i = i + 1;
end
if (Erro > toler)
    erflag=1;
else
    erflag=0;
end
modosFLAMB = v/aux;
lambFLAMB = lamb2;