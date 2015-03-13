function [R, c, dR]=con30mo(x)

R(1)=cos(x(1))-x(2);
R(2)=-x(1)+cos(x(2));
c=[];
dR=[-sin(x(1)) -1
	-1 -sin(x(2))];
