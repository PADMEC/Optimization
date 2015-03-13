function varargout = fesolrm2(munew,SN,tipo)
%
% Carry out FE analysis
%
%  load trussdata.mat
%
% First offline stage:
%
%
%  analysis
%  
% error estimator
%
 
%
% output:
%
%
%  analysis  
    global link

  area=munew(link(:,1));
  kvi=[];Z0=[];
  Lambd=[];k0n=[];
  Zd=[];kdi=[];kmni=[];
  freq=[];D=[];
 tic;
    [fn,kni,Z,kvi,Z0,kdi,kmni,Zd]  =  feoffline2(SN,tipo);
off_line_time=toc
    ti=cputime;
    [alph,u,sig,sn,vol,floc,kn,alfv,du,Lambd,k0n,D,freq] = feonline2(munew,tipo,fn,kni,Z,kvi,Z0,kdi,kmni,Zd);
    t=cputime-ti
    
    varargout={fn,kni,Z,alph,t,u,sig,sn,vol,floc,alfv,du,Lambd,D,freq};

    [co,tjjqq,vjq] = erroffline(Z);
        

%
% error estimator
%  
%   [dn]  =  erronline(alph,co,munew,ndvab,tjjqq,vjq);
   
%  plots if required
%
  iplot = 0;
  if iplot == 1 
      mun(1) = 0.
      mun(2) =  0.75
      mun(3) = 1.
      for i = 1:99
          mun(1) = mun(1) + 0.1*i;
          [snp,alphp,knp,volp] = feonline2(fn,kni,mun);
 %         
            [alphp,u,sig,snp,volp,floc] = feonline2(munew,fn,kni,Z);
 %         put everything in vectors
 %
          snv(i) = snp;
          volpv(i) = volp;
      end
      %snes = snv*10e-7
      %voles = volpv*10e-5
      %plot(snes,voles,'-')
      plot(snv,volpv,'-')
  end
  