function [e,dl,du,Aregion,ix,ipp,ipma]=UnifmComput(f,ptsp)

[npf,dim]=size(f);
mxf=max(f);
mnf=min(f);
delf=mxf-mnf;

onek=ones(npf,1);
f=(f-onek*mnf)./(onek*delf);
dl=onek*inf;
du=onek*0;
ipp=onek*0;
ipma=onek*0;
id_p=1:npf;
for ip=id_p
    jps=id_p;
    jps(ip)=[];
    for jp=jps
        %Point-Point Distance Calculation
        dist=norm(f(ip,:)-f(jp,:));
        
        if (dist<dl(ip))
            dl(ip)=dist;
            %ind d pontos proximos
            ipp(ip)=jp;
        else
        end
        fm=(f(ip,:)+f(jp,:))/2;
        distm=sum((onek*fm-f).^2,2).^.5;
        [dml]=min(distm);
        if (dml>=(dist*(.5-1e-7)))
            du(ip)=dist;
            %ind d pontos maior area
            ipma(ip)=jp;
        end
    end
end

maxd=0;
for ip=1:dim
for jp=1:dim
    dist=norm(f(ip,:)-f(jp,:));
    if (dist>maxd)
        maxd=dist;
    end
end
end
%vt=~a(maxd)^dim;
%vi=~a(dm)^dim;
%npf.vi=vt -> dm=maxd/(npf^(1/dim)) 
%4*pi*r^3/3=v
%r=(v/4/pi*3)^.5
%dm=maxd/npf;

 ix=[];
% fok=f;
% iok=id_p;
for ip=id_p
    fok=f;
    fok([ix ip],:)=[];
    dist=fok(:,1)*0;
    for id=1:dim
        dist=dist+(f(ip,id)-fok(:,id)).^2;
    end
    mdist=min(dist)^.5;
    if (mdist<1e-6)
%ind d pontos muito proximos
        ix=[ix ip];
%         iok=id_p;
%         iok(ix)=[];
    end
end

rm=maxd/2/(npf^(1/dim));
ds=[dl;du];
Ads=ds;
Ads([ix npf+ix])=0;
Ads(Ads>rm)=rm;
Aregion=sum(Ads);
e=std(ds)/mean(ds);