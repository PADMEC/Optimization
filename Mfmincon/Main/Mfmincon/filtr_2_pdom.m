function [pdom,pndom]=filtr_2_pdom(fok,FT)

[ndp,ndf]=size(fok);
[ndpT,ndfT]=size(fok);

%Normalize fok
mxf=max(fok);
mnf=min(fok);
delf=mxf-mnf;
onek=ones(ndp,1);
fok=(fok-onek*mnf)./(onek*delf);

%Normalize FT
mxf=max(FT);
mnf=min(FT);
delf=mxf-mnf;
onek=ones(ndpT,1);
FT=(FT-onek*mnf)./(onek*delf);

pdom=[];pndom=[];
npd=0;npnd=0;
for i=1:ndp
    pti=fok(i,:);
    for ifun=1:ndf
        Ri=pti(ifun);
        vecm{ifun}=find(FT(:,ifun)<(Ri-abs(Ri)/1e6));
    end
    int1=vecm{1};
    for iv=2:ndf
        int1=intersect(int1,vecm{iv});
    end
    dompnt=int1;
    if ~isempty(dompnt)
        npd=npd+1;
        pdom(npd)=i;
    else
        npnd=npnd+1;
        pndom(npnd)=i;
    end
end