function SN=amostraRB(ndv,xmin,xmax,ndb)

if nargin<4
    ndb=inputdlg('Numero de Amostragem');

    ndb=cell2mat(ndb);
    ndb=str2num(ndb);
end
ndd=ceil(log(ndb)/log(ndv));
if length(xmin)==1
    SNm=gridsamp([(xmin*4/5+xmax/5)*ones(1,ndv);(xmax4/5-xmin/5)*ones(1,ndv)*.8],ndd-1)
    SNt=gridsamp([(xmin*4/5+xmax/5)*ones(1,ndv);(xmax4/5-xmin/5)*ones(1,ndv)*.8],ndd);
elseif length(xmin)==ndv
    SNm=gridsamp([(xmin*4+xmax)/5;(xmax*4+xmin)/5],ndd-1);
    SNt=gridsamp([(xmin*4+xmax)/5;(xmax*4+xmin)/5],ndd);
end
nm=size(SNm,1);
nt=size(SNt,1);
nrest=ndb-nm;
for i=1:nm
    eq=find(sum(abs(ones(nt,1)*SNm(i,:)-SNt),2)==0);
    if length(eq)>0
        SNt(eq,:)=[];
        nt=nt-length(eq);
    end
end

%%Linearmente Dependentes
if size(SNt,1)>2*ndb*1e6
    difn=size(SNt,1)-2*ndb;
    nsn=size(SNt,1);
    for i=1:difn
        isn=ceil(rand*nsn);
        SNt(isn,:)=[];
        nsn=nsn-1;
    end
end
%[ldep,lld]=lindep(SNt);
lld=[];

if length(lld)>0
    lld(1)=[];
    nldep=length(lld);
    for i=1:nldep
        SNt(lld(i),:)=xmin+rand(1,ndv).*xmax;
    end
end

dsn=fix((ndd^ndv-1)/ndb);
nsn=size(SNt,1);
SN=SNm;
for i=1:nrest
    isn=ceil(rand*nsn);
    SN(end+1,:)=SNt(isn,:);
    SNt(isn,:)=[];
    nsn=nsn-1;
end

for i=1:ndb*ndv
    %SN(i)=SN(i)+(xmax-xmin)*rand/20;
	SN(i)=SN(i)*(1+(rand-.5)/8);
end
