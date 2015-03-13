function [sig,E,dE]=strain2stress(mat_curv,es0,E0)

if nargin==3
    E1=mat_curv(1,2)/mat_curv(1,1);
    fat=E0/E1;
end

[nelm,m] = size(es0);
M=[0 0;mat_curv;mat_curv(end,:).*[1000 10]];
ncp=size(M,1);
trac=sign(es0);
E=es0*0;
sig=E;
for icp = 2:ncp
    iel=find((abs(es0)>=M(icp-1,1))&(abs(es0)<=M(icp,1)));
    if isempty(iel)
        %i0=ncp;
        %E(iel)=0;
        %sig(iel)=M(ncp,2);
    else
        ips=[icp-1 icp];
        de=abs(es0(iel))-M(icp-1,1);
        E(iel)=diff(M(ips,2))/diff(M(ips,1));
        sig(iel)=M(icp-1,2)+E(iel).*de;
    end
    
    if nargout>2
        if i0==ncp
            E2=0;
            dE(iel)=(E2-E(iel))/diff(M(ips,1));
        else
            E2=diff(M(ips+1,2))/diff(M(ips+1,1));
            dE(iel)=(E2-E(iel))/(diff(M([i0-1 i0+1],1))/2);
        end
    end
end
sig((abs(es0)>M(ncp,1)))=M(ncp,2);
sig=sig.*trac;

if nargin==3
    sig=sig.*fat;
    E=E.*fat;
    if nargout>2
        dE=dE.*fat;
    end
    
end
    