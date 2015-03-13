function [sig,V,D]=KoskTruss(x)

%[nl,nc]=size(x);

%MN,m
F=0.020;
E=2*1e5;
L=1;
Ls=L*[2^.5 1 2];

%tt=pi/6;[cos(tt)^2 cos(tt)*sin(tt) sin(tt)^2];

%Ke=[c² sc;sc s²]
 R1=1/2*[1 -1;-1 1]/Ls(1);
 R2=1*[0 0;0 1]/Ls(2);
 R3=1/4*[3 3^.5;3^.5 1]/Ls(3);

% de cm² p/ m²
 a1=x(:,1)/1e4;a2=x(:,2)/1e4;a3=x(:,3)/1e4;
 
 %K=simple(a1*R1+a2*R2+a3*R3);
%dK=simple(det(K)*8)
%dK =a1*2^.5*a3+2*a2*a1*2^.5+a2*a3-1/2*a1*a3*6^.5;

%K1r=simple((K^-1)*dK)
% K1r =[ 2*a1+1/2*2^.5*a3,     -2*a1-1/2*2^.5*a3*3^.5
%        -2*a1-1/2*2^.5*a3*3^.5, 2*a1+4*2^.5*a2+3/2*2^.5*a3];

%E*K*d=F
%E/F*K*d=fu
%drs=simple(K1r*[1 1]')
%drs =[  a3*(1-3^.5); -a3*3^.5+8*a2+3*a3];
%d=drs/dK;
K=a1*R1+a2*R2+a3*R3;
d=K\[1 1]';
D=d*F/E;

%A=x/1e4;
%N=EA/L*R*D
%N=EA/L*R*d*F/E
%N=FA/L*R*d
%sig=N/a
%sig=F/L*R*d
Ks=F./L;

%Na=Ks(1)*R1*d/Ls(1);
%N1=norm(Na);
%c=cos(pi/4);s=sin(pi/4);
%[-c*s s^2]=[1/2 1/2]
s1=Ks*[-1 1]*d/2;
s2=Ks*[0 1]*d;
s3=Ks*[3^.5 1]*d/4;
sig=[s1,s2,s3];

V=Ls*100*x';
