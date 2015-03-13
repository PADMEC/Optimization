function [E,dE]=curv_diff(mat_curv,es0)


[m,nelm] = size(es0);
M=[0 0;mat_curv];
ncp=size(M,1);
for iel = 1:nelm
    i0=find((abs(es0(iel))<M(:,1)),1);
    if isempty(i0)
        i0=ncp;
    end
    ips=[i0-1 i0];
    E(iel)=diff(M(ips,2))/diff(M(ips,1));

    if i0==ncp
        E2=0;
        dE(iel)=(E2-E(iel))/diff(M(ips,1));
    else
        E2=diff(M(ips+1,2))/diff(M(ips+1,1));
        dE(iel)=(E2-E(iel))/(diff(M([i0-1 i0+1],2))/2);
    end
end