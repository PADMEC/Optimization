function [sig,V,d]=KoskTruss4b(x,e4)
%Messac
%[nl,nc]=size(x);


%MN,m
%KN,cm
F=0.010;
F=10;
E=2*1e5;
E=2e4*ones(1,4);
if nargin==2
    E(4)=e4;
end
L=2;
L=100;
% de cm² p/ m²
    %a=x/1e4;
    a=x;


Ls=L*[2 sqrt(2) sqrt(2) 1];
tt=[180 (180+45) -45 0];
nb=4;
vF=[0 -1 0 -1];
glb=[1 0;1 2;0 2;0 2];

%tt=pi/6;[cos(tt)^2 cos(tt)*sin(tt) sin(tt)^2];

%Ke=A/L[c² sc;-sc s²]
for i=1:nb
    c=cosd(tt(i));s=sind(tt(i));
    R{i}=[c^2 s*c;s*c s^2]*E(i)*a(i)/Ls(i);
end
%R{1}([3:4])
% 
%  R1=1/2*[1 0;0 0]/Ls(1);
%  R2=1*[0 0;0 1]/Ls(2);
%  R3=1/4*[3 3^.5;3^.5 1]/Ls(3);

%glb
    K=[R{1}+R{2} -R{2};
       -R{2} R{2}+R{3}+R{4}];

% for i=1:nb
%  K=K+a(i)*R{1};
% end

d=K\vF'*F;

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
% K=a1*R1+a2*R2+a3*R3;
% d=K\[1 1]';
% D=d*F/E;

%A=x/1e4;
%N=EA/L*R*D
%N=EA/L*R*d*F/E
%N=FA/L*R*d
%sig=N/a
%sig=F/L*R*d
%Ks=F./L;

%Na=Ks(1)*R1*d/Ls(1);
%N1=norm(Na);
%c=cos(pi/4);s=sin(pi/4);
%[-c*s s^2]=[1/2 1/2]
sig=[0,0,0,0];
for i=1:nb
    c=cosd(tt(i));s=sind(tt(i));
    if glb(i,1)==1
        sig(i)=-[c s]*d(1:2)*E(i)/Ls(i);
    else
        sig(i)=0;
    end
    if glb(i,2)==2
        sig(i)=sig(i)+[c s]*d(3:4)*E(i)/Ls(i);
    end
    
    ae = [-c -s c s];
    id = (glb(i,:)~=0);
    glbi = glb(i,id);
    gb=[glbi*2-1 glbi*2];
    uo = d(gb);
    du(i) = ae(gb)*uo;
end
k=E./Ls;
sig = du.*k;
% s1=Ks*[-1 1]*d/2;
% s2=Ks*[0 1]*d;
% s3=Ks*[3^.5 1]*d/4;
% sig=[s1,s2,s3];

V=Ls*a';
