function [T,mf,l,n,t]=MO_precomputs(ff,metodo)

[nfun,npts]=size(ff);
mxf=max(ff,[],1);
mf=min(ff,[],1);
l=mxf-mf;
l(l==0)=1;
if metodo==3||(metodo==5)
    %NBI
        
    T=ff';
    %for ifun=1:nf
    T=T-mf'*ones(1,npts);

    t=0;
    n=-sum(T,2)'/100; % t(%), -100<t<100, modified 01:2015
    %if norm(n)>0., n=n/norm(n); end
else
    %NNC
    %Normalizar
    fns=zeros(npts,nfun);
    for i=1:npts
        fns(i,:)=(ff(i,:)-mf)./l;
    end
    T=fns';
    n=zeros(npts-1,nfun);
    for i=1:npts-1
        n(i,:)=fns(end,:)-fns(i,:);
    end

    t=[];
end