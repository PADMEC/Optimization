 function kge = kgeom(esf, comp, ang, glb)
 % Monta a matriz geometrica global

 % if ndin == 3
%     Truss3D=1;
% else
%     Truss3D=0;
% end
% if Truss3D
% 	ke=[-1  0  0  1	 0  0
% 	     0 -1  0  0	 1  0
% 	     0  0 -1  0	 0  1
% 	     1  0  0 -1	 0  0
% 	     0  1  0  0	-1  0
% 	     0  0  1  0  0 -1 ];
% 
% else
% 	ke=[-1  0  1  0
% 	     0 -1  0  1
% 	     1  0 -1  0
% 	     0  1  0 -1 ];
% end

[ndin,nelm] = size(ang);
nglb = max(max(glb));
kge = sparse(nglb,nglb);

kuni=eye(ndin);
ke=[kuni  -kuni
     -kuni kuni];


% Assemble the global geometric matrix
for i = 1:nelm
    id = find(glb(:,i));
    gb = glb(id,i);
    kge(gb,gb) = kge(gb,gb) + esf(i)/comp(i)*ke(id,id);
end
kge = sparse(kge);