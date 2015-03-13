function R=con3KoskTruss(x)

%[nl,nc]=size(x);
[sig,V,D]=KoskTruss(x);

R=[];
R=[R  x-2     .1-x];
R=[R  [sig-200 -200-sig]/200];

