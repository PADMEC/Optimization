%%%%
% INTERACTIVE Differenciation Test

%K=2+x^2/2+s^2/10;
%F=20+s;
%s=1/d;

x=1;
%syms x s d real
%K=2+x^2/2+1/d^2/10;
%F=20+1/d;

%ds=solve(K*d-F,d);
ds=(5^(1/2)*(9*x^2 + 2036)^(1/2) + 100)/(5*x^2 + 20);
%dds=diff(ds(1),x);
dds=(9*5^(1/2)*x)/((5*x^2 + 20)*(9*x^2 + 2036)^(1/2)) - ...
    (10*x*(5^(1/2)*(9*x^2 + 2036)^(1/2) + 100))/(5*x^2 + 20)^2;
K=2+x^2/2+1/ds^2/2;
F=20+1/ds;
s=1/ds;

%Kd=F
%Kd'=F'-K'd
%K'= dK/dx + dK/ds*ds/dx
%K'1= dK/dx
%F'1= dF/dx
dkdx1=x;
dfdx1=0;
%Kd'=F'1-K'1d
dddx1=dfdx1-dkdx1*ds
%ds/dx
dsdx=-ds^-2;
dkds=s;

dkdx1=x+dkds*dsdx;
 