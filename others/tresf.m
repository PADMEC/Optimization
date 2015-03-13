	 function [esf,en] = tresf(area, comp, els, ang, glb, u)
% Calcula os Esforcos
     
[m,nelm] = size(comp);
esf = zeros(1,nelm);
en=esf;
for i = 1:nelm

    ae = [-cos(ang(i)) -sin(ang(i)) cos(ang(i)) sin(ang(i))];
    id = find(glb(:,i));
    gb = glb(id,i);
    uo = u(gb);
    du(i) = ae(id)*uo;
end
en=du./comp;
kl = area.*els;
esf = kl.*en;
