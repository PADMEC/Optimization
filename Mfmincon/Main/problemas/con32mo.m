function R=con32mo(x)

R(1)=x(2)^2/2+x(1)^2/2-x(3);
if R==0
    %disp('ativa')
end
R(3:5)=-x;
R(6:8)=x-50;
%R(9)=x(1)-x(2);