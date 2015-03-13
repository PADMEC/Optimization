  function [fob,gob,fre,gre]=optsolrm2(tobj,tres,xvalu,ndvab,varargin)
  %varargin{:}=fn,kni,Z,kn0,kvi
  
    global dvol lpdva
    global fs itp
    global modos

  
  %load trussdata.mat
%
% gives the values for objective function, constraints and its
% derivatives
%
%  update the primary design variable vector
%
% 
fmof=0;
if tobj == 0
    global mo
    fmof=mo(1,3);
    if fs==3|fs==5
        %%%   NBI   %%%
        global t
        t=xvalu(end);
        xvalu(end)=[];
    end
end

munew = xvalu;
	
%
% FE analysis to calculate the functions 
%

%
% 
% Sensitivity analysis for gradient evaluation
   global link glb ang comp els ro 

    props = [ang;comp;els;ro];
    %par = [ 1e-6 ; ndvab];

%
% Save all functions and their gradients
%
nmdfl=2;
    lamb=zeros(nmdfl,1);
    D=[];
    freq=[];

    [alph,u,sig,sn,vol,floc,kn,alfv,V,lamb,k0n] = feonline2(munew,modos,varargin{:});


    varargin={varargin{:},k0n,alfv,lamb};
   
    [dsn,dsig,du,dfloc,dlamb]=seonline(munew,kn,props,glb,link,modos,sig,alph,varargin{:});

%if (sum(itp~=2)==tres)&(tobj~=3)&(fmof~=1)
    %[u,sig,esf,vol,floc]=fesol(area,props,F,glb);
%    [alph,u,sig,en,vol,floc] = feonline2(munew,varargin{:});
    
    %[du,dsig,dvol,dfloc]=sensol(area,props,par,F,glb,link,lpdva,sig,u,vol,floc,gkx);
%    [dsn,dsig,du,dfloc]=seonline(munew,props,glb,link,sig,alph,varargin{:});

%    dlamb=zeros(ndvab,nmdfl);

%else
%end
%disp(floc(1:10))
    [fre,gre,fob,gob]=totim(tres,tobj,vol,dvol,sn,dsn,sig',dsig,u,du,floc',dfloc,ndvab,munew,props,lamb(1:nmdfl)',dlamb(:,1:nmdfl));    
