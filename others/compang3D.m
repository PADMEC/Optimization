% sub-rotina para calcular o comprimento e os angulos
% dos elementos.

function [ang,L] = compang3D(X,Y,Z,conect)

[nelm,n] = size(conect);
ang=zeros(3,nelm);
for i = 1:nelm
	 noj = conect(i,1);
	 nok = conect(i,2);
	 L(i) = sqrt((X(noj) - X(nok))^2 + (Y(noj) - Y(nok))^2 +...
         (Z(noj) - Z(nok))^2);

    V(1) = X(conect(i,2)) - X(conect(i,1));
    V(2) = Y(conect(i,2)) - Y(conect(i,1));
    V(3) = Z(conect(i,2)) - Z(conect(i,1));
    %ang(1,i) = atan(B(i)/A(i));
    ang(:,i)=V'/L(i);

end;

