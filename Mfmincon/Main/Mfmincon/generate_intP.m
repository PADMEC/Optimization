function [Pi]=generate_intP(P,npi)

[np,nf,nc]=size(P);
%pt=P*0;
Ps=[];
for j=1:nc
    %pt(:,:,i)=pt(:,:,i)+P(:,:,j)*weigs(j);
    Ps=[Ps;P(:,:,j)];
end
logplot=0;
if logplot==1
    figure, hold on,title('Projected Points (fix & initial)')
    if nf==2
        plot(Ps(:,1),Ps(:,2),'*')
    elseif nf>2
        scatter3(Ps(:,1),Ps(:,2),Ps(:,3),35,Ps(:,nf),'filled')
    end
end
 rest=nc*np;
%  ind=1:rest;
% 
% for j=ind
%     pis=ones(rest,1)*Ps(j,:);
%     dists(:,j)=sum((pis-Ps(ind,:)).^2,2);
% end
% 
% indd=1;
% for j=ind
%     iok(j)=ind(indd);
%     ind(ind==iok(j))=[];
%     rest=rest-1;
%     pis=ones(rest,1)*Ps(iok(j),:);
%     dists=sum((pis-Ps(ind,:)).^2,2);
%     [dmin,indd]=min(dists);
% end

if logplot==1
%[~,~,~,~,~,ipp,ipma]=UnifmComput(P,0);
 ix=[];
% fok=f;
% iok=id_p;
id_p=1:np;
ipps=zeros(np,nf-1);
for idim=1:nf-1
    for ip=id_p
        if idim>1
            ix=ipps(ip,1:(idim-1));
        else
            ix=[];
        end
        fok=P;
        idok=id_p;
        fok([ix ip],:)=[];
        idok([ix ip])=[];
        dist=fok(:,1)*0;
        for id=1:nf
            dist=dist+(P(ip,id)-fok(:,id)).^2;
        end
        [mdist,imin]=min(dist);
        %if (mdist<1e-6)
    %ind d pontos muito proximos
            ipps(ip,idim)=idok(imin);
    %         iok=id_p;
    %         iok(ix)=[];
        %end
    end
end

for ip=id_p

    if nf==2
        vx=[P(ipps(ip,:),1),ones(nf-1,1)*P(ip,1)];
        vy=[P(ipps(ip,:),2),ones(nf-1,1)*P(ip,2)];
        plot(vx,vy)
    elseif nf>2
        vx=[P(ipps(ip,:),1),ones(nf-1,1)*P(ip,1)];
        vy=[P(ipps(ip,:),2),ones(nf-1,1)*P(ip,2)];
        vz=[P(ipps(ip,:),3),ones(nf-1,1)*P(ip,3)];
        plot3(vx,vy,vz)
    end
    %plot(P(i2pp,:),P(ip,:),'r')
end
end

%y1=a+bx1
%min((x0-x1)^2+(y0-y1)^2)
%min((x0-x1)^2+(y0-a-bx1)^2)
%min((x0^2-2x0.x1 +x1^2 + y0^2-2a.y0+2a.b.x1-2y0.bx1+a^2+b^2.x1^2
%min(((-2x0+2a.b-2y0.b)x1  +(1+b^2)x1^2 + y0^2-2a.y0+a^2 +x0^2
%dR/dx1=(-x0+a.b-y0.b)2 +(1+b^2)x1.2=0
%x1=(x0-(a-y0)b)/(1+b^2)

%tp=Px-p0)*(

%min((x0^2-2x0.x1 + y0^2-2y0.(a+bx1)+a^2 + b^2.x1^2+2a.b.x1+x1^2 
%dR/dx0=2x0-2x1=0; %dR/dy0=2y0-2(a+b.x1)=0;
%x0=x1; %y0=2(a+b.x1)=0;


P0=mean(Ps);
sample_num_cvt=rest*nf*10;
[P0,cvt_time1]=uniformCVT_IntPoint(P0',sample_num_cvt,[],P',[]);

nrtest=1;

Pis=zeros(npi,nf,nrtest);
evs=zeros(1,nrtest);
for ipi=1:nrtest
    %rand0=rand(npi,2);
    rand0=lhsdesign(npi,2);
    ind=nf+ceil(rand0(:,1)*(rest-nf));
    w1=rand0(:,2)*.7;w2=1-w1;
    Pi=kron(P0,w2)+Ps(ind,:).*(w1*ones(1,nf));
    Pis(:,:,ipi)=Pi;
    %[ev]=UnifmComput(Pi,10);
    %evs(ipi)=ev;
end
bestPir=1;
%[minev,bestPir]=min(evs);
ipi=bestPir;
    Pi=Pis(:,:,ipi)';

if logplot==1
    plotnum(Pi')
end
% 
% figure, plotnum(Ps(iok,:))
% tot_weig=(nc+1)*nc/2;
% for i=1:nc
% weigs=(mod((0:nc-1)+nc-i,nc)+1)/tot_weig;
% for j=1:nc
%     pt(:,:,i)=pt(:,:,i)+P(:,:,j)*weigs(j);
% end
% plotnum(pt(:,:,i))
% end