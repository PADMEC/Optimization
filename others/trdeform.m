	 function Tu_strain = trdeform( comp, ang, glb, u)
% comput Matriz Transformation - displacement x strain
     
ngl=max(max(glb));
[m,nelm] = size(comp);
%du=zeros(nelm,1);
Tu_strain=zeros(nelm,ngl);

for i = 1:nelm

    ae = [-cos(ang(i)) -sin(ang(i)) cos(ang(i)) sin(ang(i))];
    id = find(glb(:,i));
    gb = glb(id,i);
    %uo = u(gb);
    %du(i) = ae(id)*uo;
    Tu_strain(i,gb)=Med(i,gb)+ae(id)/comp(i);
end
%en=du./comp';
%en2=Med*u;