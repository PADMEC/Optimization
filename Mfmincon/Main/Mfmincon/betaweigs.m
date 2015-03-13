function beta=betaweigs(dim,npts,mth)
%BETA - parameter create
%
% input:
% dim = output dimension
% npts = number of points
%
% output:
% beta = [npts x dim]

    %b=gridsamp([ones(1,dim-1)*0;ones(1,dim-1)*1],npts);
    b=(fullfact(ones(1,dim-1)*npts)-1)/(npts-1);
    sumb=sum(b,2);
    if mth==1
        %Include bondary
        elimb=find((sumb>1|sumb==0));
        for i=1:dim-1, elimb=[elimb;find(b(:,i)==1)];end
    else
        %only int points
        elimb=find((sumb>=1|sumb==0));
        for i=1:dim-1, elimb=[elimb;find(b(:,i)==1|b(:,i)==0)];end
    end
    b(elimb,:)=[];
    sumb(elimb)=[];
    br=1-sumb;
    beta=[b br];
    
end
%BETA - parameter create
%b2=gridsamp([zeros(1,paret_surf_dim-1);ones(1,paret_surf_dim-1)],ndiv);
%nd=ptsp;
%ibs=0:1/nd:1;
%b2=ibs';
%sumb=sum(b2,2);
%elimb=find((sumb>1|sumb==0));
%for i=1:probdim(1)-1, elimb=[elimb;find(b2(:,i)==1)];end
%b2(elimb,:)=[];
%sumb(elimb)=[];
%br=1-sumb;
%b2=[b2 br];
%%%%%%%%%%%%%%%%%%%%