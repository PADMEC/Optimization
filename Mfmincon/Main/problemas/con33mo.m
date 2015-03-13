function R=con33mo(x)

R(:,1)=-5+x(:,2).^2+x(:,1).^2-x(:,3);
if R==0
    %disp('ativa')
end
%a=2500;
%c=100;
%R(2)=-1;%a+c/(c/a+x(1)+x(2));
%if R(2)==0
%    disp('ativa')
%end
R(:,2:4)=-x;
R(:,5)=x(:,3)-(x(:,1)+x(:,2))*5;
%R(7)=-x(1)-x(2)+7;
%R(7)=x(1)-x(2);