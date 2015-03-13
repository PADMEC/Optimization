 function kgb = gkeglb(area, comp, els, ang, glb,idvab,lpdva,link)  
 % Monta a matriz de rigidez global

[m,nelm] = size(comp);
nglb = max(max(glb));
kgb = zeros(nglb,nglb);
ielem = lpdva(idvab);

for jelem = 1:nelm
   kdvab = link(jelem,1); 
   kelem = lpdva(kdvab);
    
   ke = gkeltr(area(jelem),comp(jelem),els(jelem),ang(jelem),kelem,jelem,ielem);
   id = find(glb(:,jelem));
   gb = glb(id,jelem);
   kgb(gb,gb) = kgb(gb,gb) + ke(id,id);
end;
