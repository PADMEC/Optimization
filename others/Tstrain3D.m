	 function Mu_strain = Tstrain3D( comp, cosn, glb)
% comput Matriz Transformation - displacement x strain
% strain = M*u;
     
ngl=max(max(glb));
nelm = length(comp);
%du=zeros(nelm,1);
Mu_strain=zeros(nelm,ngl);

for i = 1:nelm

	ae = [-cosn(1,i)  -cosn(2,i)  -cosn(3,i)  cosn(1,i)  cosn(2,i)  cosn(3,i)];
    id = find(glb(:,i));
    gb = glb(id,i);
    Mu_strain(i,gb)=Mu_strain(i,gb)+ae(id)/comp(i);
    %uo = u(gb);
    %du(i) = ae(id)*uo;
end

%en=du./comp';
%en2=Med*u;