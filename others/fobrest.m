function [fde,fte,tpres]=fobrest(arq,ngl,area,els,ro,comp,cm,em,d)

	global con id idm kb MinLambd

con=[];
fde=0;fte=0;
d=[];kb=[];idm=[];
if isfield(arq,'d')
    if length(arq.d)==1
        d=arq.d*ones(1,ngl);
        id=find(arq.d);
    else
        id=find(arq.d);
        d=arq.d(id);
    end
end
if isfield(arq,'idm')
    fde=length(arq.idm);
end
if isfield(arq,'kb')
    kb=arq.kb;
    fte=length(arq.kb);
end

if isfield(arq,'idm')
    idm=arq.idm;
end
if isfield(arq,'cm')
    cm=arq.cm;
end
if isfield(arq,'em')
    em=arq.em;
end

global itp SIGu
itp=find(arq.tpres);
nel=length(area);
clb=[];cub=[];
for i=1:sum(arq.tpres)
    if itp(i)==1 %Volume
        if isfield(arq,'vol0')
            clb=arq.vol0*0;
            cub=arq.vol0;
        else
            clb=sum(area.*comp.*ro)*.99999;%*ro(1);%   VOL
            cub=sum(area.*comp.*ro);%*ro(1);% VOL
        end

   elseif itp(i)==5 %Flambagem Global
       %instability!
       
       veclamb=ones(1,min(5,ngl));
       if isempty(MinLambd)
           clb=1*veclamb;%ceil(nglb/4));
       else
           clb=MinLambd*veclamb;
       end
       cub=inf*veclamb;%ceil(nglb/4));

    elseif itp(i)==0 % Flambagem Local
        E=els;
        clb=E-E-inf;%-E*pi/4./comp.^2;
        cub=E-E+1;%1e9*E*pi/4./comp.^2;

    elseif itp(i)==3
        if isfield(arq,'cm')
            clb=arq.cm*ones(1,nel);
        end
        if isfield(arq,'em')
            cub=arq.em*ones(1,nel);
            disp('[clb,cub] = arq.cm,em')
            disp([arq.cm, arq.em])
        else
            clb=-SIGu/.99999*ones(1,nel);
            cub=SIGu*.99999*ones(1,nel);
        end

    elseif itp(i)==2
        disp('[-clb,cub] = arq.d')
        disp([arq.d])

        clb=-d;
        cub=d;
        
    end
    con = [con [clb;cub]];
end

tpres=(itp);
if isempty(itp), tpres=0; end
