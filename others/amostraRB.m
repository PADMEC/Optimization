function [SN,SNm]=amostraRB(ndv,xmin,xmax,ndb)

if nargin<4
    ndb=inputdlg('Numero de Amostragem');
    ndb=cell2mat(ndb);
    ndb=str2num(ndb);
    if nargin==1
        arq=ndv;
        ndv=max(arq.arex);
        xmin=arq.xmin(1);
        xmax=arq.xmax(1);
    end
end
ndd=ceil(exp(log(ndb)/ndv));
dv=(xmax-xmin)/200;
%dv=(xmax-xmax);
if length(xmin)==1
    SNm=gridsamp([(xmin+dv)*ones(1,ndv);(xmax-dv)*ones(1,ndv)],ndd-1);
    SNt=gridsamp([(xmin+dv)*ones(1,ndv);(xmax-dv)*ones(1,ndv)],ndd);
    xmax=xmax*ones(1,ndv);
    xmin=xmin*ones(1,ndv);
else
    SNm=gridsamp([(xmin+dv);(xmax-dv)],ndd-1);
    SNt=gridsamp([(xmin+dv);(xmax-dv)],ndd+2);
end
if ndb==ndd^ndv
    SNm=gridsamp([(xmin+dv);(xmax-dv)],ndd);
end
SN=SNm;
nm=size(SNm,1);
nt=size(SNt,1);
nrest=ndb-nm;

ond=ones(nt,1); distNm=[];
for i=1:nm
    distNm(:,i)=sum(((ond*SNm(i,:)-SNt)./(ond*(xmax-xmin))).^2,2);
end
[mdm,mlin]=min(distNm,[],2);
[mm,mcol]=sort(-mdm);
SNt=SNt(mcol(1:nrest),:);

%%Linearmente Dependentes
%[ldep,lld]=lindep(SNt);
lld=[];
if length(lld)>0
    lld(1)=[];
    nldep=length(lld);
    for i=1:nldep
        SNt(lld(i),:)=xmin+rand(1,ndv).*xmax;
    end
end

%dsn=fix((ndd^ndv-1)/ndb);
%nsn=ndd^ndv;
% for i=1:nrest
%     isn=ceil(rand*nt);
%     SN(end+1,:)=SNt(isn,:);
%     SNt(isn,:)=[];
%     nt=nt-1;
% end
%
% dn=floor(nt/nrest);
% nt/6
% nts=floor([nt/6 nt/3 nt/2])+floor(mod(nt,[6 3 2])/3);
% dns=round(nts/nrest*3);
% rest=mod(nts,nrest/3);
% i1=[              1:dns(1):sum(nts(1:1))-rest(1)];
% i2=[       nts(1)+1:dns(2):sum(nts(1:2))-rest(2)];
% i3=[sum(nts(1:2))+1:dns(3):sum(nts(1:3))-rest(3)];
%SN=[SN;SNt([i1 i2 i3],:)];
% nti=floor(nt/3);
% nreti=floor(nrest/3)*ones(3,1);
% rti=mod(nrest,3);
% nreti(1:rti)=nreti(1:rti)+1;
% for i=1:3
%     
%     SN=[SN;SNt(1:nreti(i),:)];
%     SNt(1:nti,:)=[];
% end

%an+bn^2=nt
%a+b=1
if nrest==1
    in=1;
elseif nrest>1
    ab=[nrest nrest^2;1 1]\[nt 1]';
    i2=ab(1)*2+ab(2)*4;
    if i2>1.5
        in=round(ab(1)*(1:nrest)+ab(2)*(1:nrest).^2);
    else
        in=ceil(ab(1)*(1:nrest)+ab(2)*(1:nrest).^2);
    end
end
%in(end)=nt;
% prox=0;
% while length(prox)>0
%     prox=find(diff(in)==0);
%     in(prox)=in(prox)+1;
% end
SN=[SN;SNt(1:nrest,:)];
nsn=size(SN,1);
nrest=ndb-nsn;
SN=[SN;SNt(1:nrest,:)];
for i=1:ndb
    for j=(i+1):ndb
        if norm(SN(i,:))==0
            dd=0;
        else
            dd=SN(i,:)'\SN(j,:)';
        end
        if norm(SN(j,:))==0, errl=norm((SN(i,:)*dd-SN(j,:)));else
        errl=norm((SN(i,:)*dd-SN(j,:)))/norm(SN(j,:));
        end
        if errl<1e-6
            %rdi=rand(1,ndb)/10;
            rdi=[1:ndv]/10/ndv;
            SN(i,:)=SN(i,:).*(1-rdi)+rdi.*(xmax+xmin)/2;
            %%%%%%%%%%%PARA converg rate test%%%%%%%%%%%%%%%%%%%%%%%%%%%
            SN(j,:)=SN(j,:).*(1-rdi)+rdi.*(xmax+xmin)/2;
        end
    end
end
%SN=[SN;SNt(1:dn:nt-rest,:)];
% for i=1:ndb*ndv
%     %SN(i)=SN(i)+(xmin/xmax)*rand/2;
% 	SN(i)=SN(i)*(1+(.5-rand)/10);
% end