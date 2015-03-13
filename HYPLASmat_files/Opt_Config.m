function [fde,fte,tpres]=Opt_Config(tpres, topt,varargin)
% based on [Fs,Us,Sigs,ien,vol  ]
         global con
         global  id
         global idm
         global kb
        global itp SIGu
xxxxxxxxxxxxxxxxxxxxx
con=[];
fde=0;fte=0;
d=[];kb=[];idm=[];


itp=find(tpres);
nel=nn;
clb=[];cub=[];
for i=1:sum(tpres)
    
    [mx,km]=max(abs(varargin{i}));
    
    
    if itp(i)==1 %Volume
        if isfield(arq,'vol0')
            clb=arq.vol0*0;
            cub=arq.vol0;
        else
            clb=sum(area.*comp.*ro)*0;%*ro(1);%   VOL
            cub=sum(area.*comp.*ro);%*ro(1);% VOL
        end

   elseif itp(i)==-1 %Flambagem Global
        clb=ones(1,2);%ceil(nglb/4));
        cub=inf*ones(1,2);%ceil(nglb/4));

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
        clb=-d;
        cub=d;
        
    end
    con = [con [clb;cub]];
end

tpres=(itp);
if isempty(itp), tpres=0; end
