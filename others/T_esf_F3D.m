	 function Mesf_F = T_esf_F3D( cosn, glb)
% comput Matriz Transformation - esf -> F
% F = M*esf;

ngl=max(max(glb));
nelm = size(glb,2);
Mesf_F=zeros(ngl,nelm);

for i = 1:nelm
    ae = [-cosn(:,i); cosn(:,i)];
    id = find(glb(:,i));
    gb = glb(id,i);
    Mesf_F(gb,i)=Mesf_F(gb,i)+ae(id);
    
    %Fi = ae(id)*esf(i);
    %F(gb)=F(gb)+Fi;
    %kl = area(i)*els(i)/comp(i);
end

%F2=Mesf_F*esf;
