function [X,cpu_time]=uniformCVT_IntPoint(X,sample_num_cvt,Vr,xvr,P0)
%
%X=x0
%p=dimension
%n=size x
%sample_num_cvt=auxiliar points number
%
global toll
if (isempty(toll))
    toll=1e-6;
elseif (toll==0)
    toll=1e-6;
end
aaa=size(xvr);
p=aaa(1);
n=size(X,2);
% for i=1:aaa(3)
%     Xfix=[Xfix; fs(:,:,i)];
%     mns(i,:)=min(fs(:,:,i));
%     mxs(i,:)=max(fs(:,:,i));
% end

if length(aaa)<3
    i=1;
    Xfix=xvr(:,:,i)';
    mns(i,:)=min(xvr(:,:,i),[],2)';
    mxs(i,:)=max(xvr(:,:,i),[],2)';
    nfix=aaa(2);
else
    Xfix=[];
    for i=1:aaa(3)
        Xfix=[Xfix; xvr(:,:,i)'];
        mns(i,:)=min(xvr(:,:,i),[],2)';
        mxs(i,:)=max(xvr(:,:,i),[],2)';
    end
    nfix=prod(aaa(2:3));
end

n=n+nfix;
mns=min(mns,[],1);
mxs=max(mxs,[],1);
X=[Xfix; X'];

%number of divisions in uniform grid points
ndd=round(exp(log(sample_num_cvt)/p));
%ndd^p=sample_num_cvt
nddm=ceil(exp(log(n)/p));
if nddm>ndd
    ndd=nddm;
end
auxP=gridsamp([mns;mxs],ndd);
% 
logplot=0;
%voronoi(X(:,1),X(:,2),X(:,3));

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
nt=size(auxP,1);

Kfix=convhulln([X]);
out_auxP=zeros(nt,1);
for i=1:nt
    Khull=convhulln([auxP(i,:);X(Kfix(:,1),:)]);
    if ~isempty(find(Khull==1,1))
        out_auxP(i)=1;
    end
end
% ix= out_auxP;
auxP(out_auxP==1,:)=[];


if logplot==1
if p==2
    figure,hold on
    plot(X(:,1),X(:,2),'.');
    plot(auxP(:,1),auxP(:,2),'og');
elseif p==3
    figure,hold on
    plot3(X(:,1),X(:,2),X(:,3),'.');
    plot3(auxP(:,1),auxP(:,2),auxP(:,3),'og');
end
end

nt=size(auxP,1);
ond=ones(nt,1); distNm=[];% distNm(i,j) = dist( auxP(i) to X(j) )
tc0=cputime;
%onefix=ones(nfix,1);
dauxP=zeros(1,500);
distNm=zeros(nt,p);
% delx=mxs-mns;
for iter=1:1000
    %auxP=lhsdesign(nt*10,p);
    %auxP=rand(nt*10,p);
    %auxP=ond*mns+auxP.*(ond*delx);
    for i=1:n
        distNm(:,i)=sum((ond*X(i,:)-auxP).^2,2);
    end
    
    [mdm,mlin]=min(distNm,[],2);
    newX=X;
    for i=(nfix+1):n
        %points of Xi region
        pvr=find(mlin==i);
        npvr=length(pvr);
        M(i)=npvr/nt-1/n;
        if npvr~=0
            dx=sum(auxP(pvr,:),1);
            %Gx(i,:)=dx/npvr;
            
            newX(i,:)=(X(i,:)+dx)/(npvr+1);
            %newX(i,:)=dx/(npvr);
            
            %%Test for point out of pareto domain
            %auxV=onefix*newX(i,:)-Xfix;
            %ordering
            %[sdist,nersts]=sort(sum(auxV.^2,2));
            %auxV*Lamb=0 ->Convex
            %V*null(V)=0, 
            %nerstV=auxV(nersts(1:(p+1)),:)';
            %P=null(nerstV);
            %abs(sum(sign(P)));
            %if abs(sum(sign(P)))==p
            %    newX(i,:)=X(i,:);
            %end
        end
        
        if logplot==1
        %%%%%%%%%%%%%%%%%%%%
        %%% PLOT REGIONS %%%
        switch i
            case (nfix+1)
            colr='k';
            case (nfix+2)
            colr='m';
            case (nfix+3)
            colr='g';
            case (nfix+4)
            colr='r';
        end
        if (i<(nfix+5))            
             plot(X(i,1),X(i,2),['.' colr]);
             plot(auxP(pvr,1),auxP(pvr,2),['x' colr]);
        end
        %%%%%%%%%%%%%%%%
        end
    end
%     figure,hold on
%      plot(X(i,1),X(i,2),'.k');
%      plot(auxP(pvr,1),auxP(pvr,2),'xk');

%     plot(hm1o,M);
%     voronoi(hvoroi,X(:,1),X(:,2));
%     hold on
%     plot(newX(:,1),newX(:,2),'.k')
%     plot(Gx(:,1),Gx(:,2),'.g')
%     voromuv(intr)=getframe(voroi);
    %hold off
    %auxP=rand(nt,p);
    dauxP(iter)=norm(newX-X)/norm(newX);
    %M1oi(iter)=norm(M)/n;
    X=newX;
    
    if logplot==1
    if iter==1
        plot(X(:,1),X(:,2),'.k');
        plot(X(i,1),X(i,2),'ok');
        plot(auxP(pvr,1),auxP(pvr,2),'x');
    elseif iter==2
        figure,hold on
        plot(X(:,1),X(:,2),'.');
        plot(X(i,1),X(i,2),'ok');
        plot(auxP(pvr,1),auxP(pvr,2),'x');
    elseif mod(iter,5)==0
        plot(X(:,1),X(:,2),'.');
        plot(X(i,1),X(i,2),'ok');
        plot(auxP(pvr,1),auxP(pvr,2),'x');
    end
    end
        
    if dauxP(iter)<toll
        %if norm(M)<toll*10
        break
        %end
    end
end
tc1=cputime;
cpu_time=tc1-tc0;
%figure
%voronoi(X(:,1),X(:,2)),hold on
%plot(Gx(:,1),Gx(:,2),'.g')
fprintf ( 1, 'Total steps: %4d, Change: %10f, norm(M)%10f , CPU Time: %f\n', iter, dauxP(iter),norm(M), cpu_time);
% figure,plot(dauxP),title('step length')
% figure,plot(M1oi),title('region size')
% figure
% movie(voromuv)

if logplot==1
%         %%%%%%%%%%%%%%%%%%%%
%         %%% PLOT REGIONS %%%
for i=(nfix+1):n
        %points of Xi region
        pvr=find(mlin==i);
        %npvr=length(pvr);
        %M(i)=npvr/nt-1/n;      
        switch (mod(i,5)+1)
            case (1)
            colr='k';
            case (2)
            colr='m';
            case (3)
            colr='g';
            case (4)
            colr='r';
            case (5)
            colr='y';
        end
        %if (i<(nfix+5))            
             plot(X(i,1),X(i,2),['+' colr]);
             plot(auxP(pvr,1),auxP(pvr,2),['.' colr]);
        %end
end
end
%         %%%%%%%%%%%%%%%%
X=X((nfix+1):end,:);