function R=con3KimA(x)

[nl,nc]=size(x);
potM=ones(nl,1)*[4 3 2];
R= (x.^potM)*[1;2;5] -1;
R=[R -x];