% sub-rotina para calcular o comprimento e os angulos
% dos elementos.

function [teta,L] = compang(X,Y,conect)

[nelm,n] = size(conect);

for i = 1:nelm
    noj = conect(i,1);
    nok = conect(i,2);
    L(i) = sqrt((X(noj) - X(nok))^2 + (Y(noj) - Y(nok))^2);

    A(i) = X(conect(i,2)) - X(conect(i,1));
    if A(i)==0
        A(i)=10^(-9);
    end
    B(i) = Y(conect(i,2)) - Y(conect(i,1));
    teta(i) = atan(B(i)/A(i));
    if A(i)<0
       teta(i)=teta(i)+pi;
    end
end

