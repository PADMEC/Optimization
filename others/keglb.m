 function [kgb,igl_ok, bars_out,gbs] = keglb(area, comp, els, ang, glb,Kelem_new)  
 % Monta a matriz de rigidez global

ndin = size(ang,1);
if ndin == 3
    Truss3D=1;
else
    Truss3D=0;
end

[~,nelm] = size(comp);
nglb = max(max(glb));


%Failure Bars
% rupelem=find(els==0);
% 
% %Graus de Liberdade das barras
% glout=glb(:,rupelem);
% glout=glout(find(glout))';
% nglo=length(glout);
% 
 out_gl=zeros(1,nglb);
% glbtest=glb;
% glbtest(:,rupelem)=[];
% barstest(rupelem)=[];
 bars_out=[];
% for igl=glout
%     nbar=find(glbtest==igl,1);
%     if length(nbar)<2
%         out_gl(igl)=1;
%         ibar=ceil(nbar/4);
%         bars_out=[bars_out barstest(ibar)];
%         
%     end
% end
% 
 igl_ok=find(out_gl==0);
 
%global sigs0  ellcurv

kgb = sparse(nglb,nglb);
%kgb=zeros(nglb);
gbs=zeros(1,nglb);

% Get the element that needed update
barstest=1:nelm;
 if Kelem_new==0
 else
    barstest=barstest(Kelem_new);
 end
for i = barstest
    %if isempty(intersect(rupelem,i))
    if Truss3D
        ke = keltr3D(area(i),comp(i),els(i),ang(:,i));
    else
        ke = keltr  (area(i),comp(i),els(i),ang(:,i));
    end
        id = find(glb(:,i));
        gb = glb(id,i);
        gbs(gb)=1;
        kgb(gb,gb) = kgb(gb,gb) + ke(id,id);
    %end
end
gbs=find(gbs);
kgb = sparse(kgb);