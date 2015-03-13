	 function Mesf_F = T_esf_F( ang, glb)
% comput Matriz Transformation - esf -> F
% F = M*esf;

ngl=max(max(glb));
[ndim,nelm] = size(ang);
Mesf_F=zeros(ngl,nelm);

if ndim==1
    %2D
    cosn=[cos(ang); sin(ang)];
elseif ndim==3
    cosn=ang;
end
for i = 1:nelm

    ae = [-cosn(:,i); cosn(:,i)];
    %disp('test 2D... T_esf_F')
    %ae = [-cos(ang(i)) -sin(ang(i)) cos(ang(i)) sin(ang(i))]';
    id = find(glb(:,i));
    gb = glb(id,i);    
    Mesf_F(gb,i)=Mesf_F(gb,i)+ae(id);
    
    %Fi = ae(id)*esf(i);
    %F(gb)=F(gb)+Fi;
    %kl = area(i)*els(i)/comp(i);
end

%F2=Mesf_F*esf;
