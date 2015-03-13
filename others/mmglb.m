function mgb = mmglb(area,comp,glb,ro)  
%###############################################################
% Monta a Matriz de Massa Global
%###############################################################

[ngln,nelm] = size(glb);
nglb = max(max(glb));
mgb = sparse(nglb,nglb);

if ngln==4
    % Truss 2D
    me=[2 0 1 0
        0 2 0 1
        1 0 2 0
        0 1 0 2];
elseif ngln==6
    % Truss 3D
    me=[2 0 1 0	1 0
        0 2 0 1	0 1
        1 0 2 0 1 0
        0 1 0 2 0 1
        1 0 1 0 2 0
        0 1 0 1 0 2];
end

for i = 1:nelm
        ml = ro(i)*area(i)*comp(i)/6;
        id = find(glb(:,i));
        gb = glb(id,i);
        mgb(gb,gb) = mgb(gb,gb) + ml*me(id,id);
end
