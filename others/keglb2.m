 function [kgb,igl_ok, bars_out] = keglb(area, comp, els, ang, glb)  
 % Monta a matriz de rigidez global

[m,nelm] = size(comp);
nglb = max(max(glb));


%Failure Bars
 rupelem=(els==0);
 bars_out=[];
 if sum(rupelem)>0
% 
% %Graus de Liberdade das barras
 glout=glb(:,rupelem);
 glout=unique(glout(:))';
 %glout=glout(glout)';
 nglo=length(glout);
 
 glbtest=glb;
 glbtest(:,rupelem)=[];
 glbtest=glbtest(:);
 %intersect()
% 
 out_gl=zeros(1,nglb);
 barstest=1:nelm;
 barstest(rupelem)=[];
for igl=glout
    nbar=find(glbtest==igl,1);
    if length(nbar)<2
        out_gl(igl)=1;
        ibar=ceil(nbar/4);
        bars_out=[bars_out barstest(ibar)];
        
    end
end
 igl_ok=find(out_gl==0);
else
     igl_ok=1:nglb;
end
 
global sigs0  ellcurv

kgb = sparse(nglb,nglb);
%kgb=zeros(nglb);

for i = barstest
    %if isempty(intersect(rupelem,i))
        ke = keltr(area(i),comp(i),els(i),ang(i));
        id = find(glb(:,i));
        gb = glb(id,i);
        kgb(gb,gb) = kgb(gb,gb) + ke(id,id);
    %end
end;
kgb = sparse(kgb);