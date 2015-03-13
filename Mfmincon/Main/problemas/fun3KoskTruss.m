function f=fun3KoskTruss(x)

%[nl,nc]=size(x);
[sig,V,D]=KoskTruss(x);

%papern=1 - Kim: adapt WS
%papern=2 - Kim: ENNC
papern=1;
if papern==1
    f=[V, abs(sig(1)), abs(sig(3))];
    
else
    m1=[3 1]*D/4;

    for i=1:3
        a=x(i);
        if a<0.9
            cp=-(a-.1)^2/.128+20;
        elseif a<1.5
            cp=3.33*a^2 + 9.67;
        elseif a<3
            cp=-(a-3)^2/.18 + 28;
        end
        C(i)=cp;
    end


    m2=V*1.5*7850/100^3 + sum(C);
    m3=V/10^3;%*7850

    f=[m1, m2, m3];
end