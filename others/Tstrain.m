	 function Mu_strain = Tstrain( comp, ang, glb)
% comput Matriz Transformation - displacement x strain
% strain = M*u;
     
ngl=max(max(glb));
[m,nelm] = size(comp);
%du=zeros(nelm,1);
Mu_strain=zeros(nelm,ngl);

for i = 1:nelm

    ae = [-cos(ang(i)) -sin(ang(i)) cos(ang(i)) sin(ang(i))];
    id = find(glb(:,i));
    gb = glb(id,i);
    Mu_strain(i,gb)=Mu_strain(i,gb)+ae(id)/comp(i);
    %uo = u(gb);
    %du(i) = ae(id)*uo;
end
%en=du./comp';
%en2=Med*u;