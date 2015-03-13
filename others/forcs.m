	 function [F] = forcs(area, ang, glb, sig)
% comput Matriz Transformation - esf -> F
% F = M*esf;

[m,nelm] = size(glb);
ngl=max(max(glb));
F = zeros(ngl,1);
%en=esf;
esf = area'.*sig;

Mesf_F=zeros(ngl,nelm);
for i = 1:nelm

    ae = [-cos(ang(i)) -sin(ang(i)) cos(ang(i)) sin(ang(i))]';

    id = find(glb(:,i));
    gb = glb(id,i);
    Fi = ae(id)*esf(i);
    
    Mesf_F(gb,i)=Mesf_F(gb,i)+ae(id);
    
    F(gb)=F(gb)+Fi;
    %kl = area(i)*els(i)/comp(i);
end
