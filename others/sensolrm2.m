function [t,du,dsig,dvol,den,gfloc,lambgrad]=sensolrm2(munew,SN,modos)
%####################################################################
%# sensitivity analysis--- direct analytical procedure	    #
%####################################################################
%
%
%--------------------------------------------------------------------
%
%   escolha um dos m‚todos:
% 
%   diferencas finitas      metd=1
%   direto analitico        metd=2
%   direto semi-analitico   metd=3
%   adjunto                 metd=4
%
%--------------------------------------------------------------------
%global dvol;
%global comp1;
%global comp2;
%global comp3;
%metd=4;

%if metd==1
%	[ugra,sigra,volgra]=senfd(area,props,par,fext,glb,link,lpdva,sig,u,vol);
%	elseif metd==2
%		[ugra,sigra,volgra]=senal(area,props,par,fext,glb,link,lpdva,sig,u,vol);
%elseif metd==3
%			[ugra,sigra,volgra]=sensal(area,props,par,fext,glb,link,lpdva,sig,u,vol);
%			else
%			[ugra,sigra,volgra]=sensaj(area,props,par,fext,glb,link,lpdva,sig,u,vol);
%
%end

  load trussdata.mat
%
% First offline stage:
%
%
%  analysis
global link glb ang comp els ro 
props = [ang;comp;els;ro];
area=munew(link(:,1));

    [fn,kni,Z,kvi,Z0]  =  feoffline2(SN,modos);

    [alph,u,sig,sn,vol,floc,kn,alfv,du,Lambd,k0n] = feonline2(munew,modos,fn,kni,Z,kvi,Z0);
    
    tic;
    %[ugra,sigra,volgra,flocgrad,lambgrad]=senal(area,props,par,fext,glb,link,lpdva,sig,u,vol,floc,gkx,V,Lambd);
    [den,dsig,du,gfloc,lambgrad]=seonline(munew,kn,props,glb,link,modos,sig,alph,fn,kni,Z,kvi,Z0,k0n,alfv,Lambd);
    t=toc;


%[dsn] = seonline(alph,fn,kn,kni,ndvab);

%  
% error estimator
%  
 %  [dn]  =  erronline(alph,co,munew,ndvab,tjjqq,vjq);
   
 %          [arq.idis,arq.ites,esf,vol,floc,du,Lambd] 

