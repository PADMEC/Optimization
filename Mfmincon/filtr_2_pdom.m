function [pndom,whodom]=filtr_2_pdom(FT)

[ndp,ndf]=size(FT);

%Normalize FT
mxf=max(FT);
mnf=min(FT);
delf=mxf-mnf;
onek=ones(ndp,1);
FT=(FT-onek*mnf)./(onek*delf);
dr=1e-4;

pndom=[];whodom=[];
iT=1:ndp;
npd=1;npnd=0;npdt=0;
while npd>0
    npd=0;
    pdom=[];
    for i=1:ndp
        pti=FT(i,:);
        for ifun=1:ndf
            Ri=pti(ifun);
            vecm{ifun}=find(FT(:,ifun)<(Ri-dr));
        end
        int1=vecm{1};
        for iv=2:ndf
            if isempty(int1)
                break
            end
            int1=intersect(int1,vecm{iv});
        end
        % Tolerance in F difference in the utopic direction 2%
        maxdf=0.02;
        dompnt=0;
        for ip=int1'
            df=sum(pti-FT(ip,:))*.70710678;
            if df>maxdf
                dompnt=ip;
                maxdf=df;
            end
        end

        if dompnt
            npd=npd+1;
            pdom(npd)=i;
            npdt=npdt+1;
            whodom(npdt,:)=[iT(i),iT(dompnt)];
        else
            npnd=npnd+1;
            pndom(npnd)=iT(i);
        end
    end
    iT(pdom)=[];
    FT(pdom,:)=[];
    [ndp]=size(FT,1);
end