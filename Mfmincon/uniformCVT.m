function X=uniformCVT(X,p,n,sample_num_cvt)
%
%X=x0
%p=dimension
%n=size x
%sample_num_cvt=auxiliar points number
%
global toll

ndd=round(exp(log(sample_num_cvt)/p));
%ndd^p=sample_num_cvt
nddm=ceil(exp(log(n)/p));
if nddm>ndd
    ndd=nddm;
end
Xt=gridsamp([zeros(1,p);ones(1,p)],ndd);
% 
% figure,hold on
% plotnum(X(:,1) ,X(:,2) ,'or')
% plotnum(Xt(:,1),Xt(:,2),'.')

figure,hold on
voronoi(X(:,1),X(:,2));
%[V,C] = voronoin(X);
%for in=1:n
%    ci=C{in};
%    ci(find(ci==1))=[];
%    plot(mean(V(ci,1)),mean(V(ci,2)),'.g')
%end
%plotnum(X(:,1),X(:,2),'.')
%xlim([0 1]),ylim([0 1])

% figure,hm1o=axes;hold on
% figure,hvoroi=axes;xlim([0 1]),ylim([0 1])

nt=size(Xt,1);
ond=ones(nt,1); distNm=[];% distNm(i,j) = dist( Xt(i) to X(j) )
tic
for intr=1:100
    for i=1:n
        distNm(:,i)=sum((ond*X(i,:)-Xt).^2,2);
    end
    [mdm,mlin]=min(distNm,[],2);
    newX=X;
    for i=1:n
        pvr=find(mlin==i);
        npvr=length(pvr);
        M(i)=npvr/nt-1/n;
        if npvr~=0
            dx=sum(Xt(pvr,:));
            %Gx(i,:)=dx/npvr;
            newX(i,:)=(X(i,:)+dx)/(npvr+1);
            %newX(i,:)=dx/(npvr);
        end  
    end
%     plot(hm1o,M);
%     voronoi(hvoroi,X(:,1),X(:,2));
%     hold on
%     plot(newX(:,1),newX(:,2),'.k')
%     plot(Gx(:,1),Gx(:,2),'.g')
%     voromuv(intr)=getframe(voroi);
    %hold off
    %Xt=rand(nt,p);
    dXt(intr)=norm(newX-X)/n;
    M1oi(intr)=norm(M)/n;
    X=newX;
    if dXt(intr)<toll
        break
    end
end
toc
figure
voronoi(X(:,1),X(:,2)),hold on
%plot(Gx(:,1),Gx(:,2),'.g')
fprintf ( 1, 'Total steps: %4d, Chenge: %10f, norm(M)%10f\n', intr, dXt(intr),norm(M) );
figure,plot(dXt)
figure,plot(M1oi)

% figure
% movie(voromuv)
