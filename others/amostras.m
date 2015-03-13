function SN=amostras(arq)

ndv=max(arq.arex);
ndb=inputdlg('Numero de Amostragem');

ndb=cell2mat(ndb);
ndb=str2num(ndb);
ndd=ceil(log(ndb)/log(ndv));
SNt=gridsamp([arq.xmin*ones(1,ndv)*1.2;arq.xmax*ones(1,ndv)*.8],ndd);

%%Linearmente Dependentes
[ldep,lld]=lindep(SNt);
lld(1)=[];
nldep=length(lld);
SNt(lld,:)=arq.xmin+rand(nldep,3)*arq.xmax;


dsn=fix((ndd^ndv-1)/ndb);
nsn=ndd^ndv;
for i=1:ndb
    isn=ceil(rand*nsn);
    SN(i,:)=SNt(isn,:);
    SNt(isn,:)=[];
    nsn=nsn-1;
end

for i=1:ndb*ndv
    %SN(i)=SN(i)+(arq.xmin/arq.xmax)*rand/2;
	SN(i)=SN(i)*(1+rand/4);
end
