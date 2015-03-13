% Realisa Análise de sensibilidades usando uma das tecnicas

  function [den,ugra,sigra,volgra,flocgrad,lambgrad]=sensol(area,props,par,fext,glb,link,lpdva,sig,u,vol,floc,gkx,modos,V,Lambd)

%####################################################################
%# Realisa Análise de sensibilidades usando uma das tecnicas 	    #
%####################################################################
%
%--------------------------------------------------------------------
%
%   escolha um dos metodos:
% 
%   diferencas finitas      metd=1  -1o.
%   direto analitico        metd=2  -2o.
%   direto semi-analitico   metd=3  -3o.
%   adjunto                 metd=4  -4o.
%
%--------------------------------------------------------------------

metd=2;

if metd==1
        [den,ugra,sigra,volgra,flocgrad,lambgrad]=senfd(area,props,par,fext,glb,link,lpdva,sig,u,vol,floc,gkx,modos,V,Lambd);

elseif metd==2
        [den,ugra,sigra,volgra,flocgrad,lambgrad]=senal(area,props,par,fext,glb,link,lpdva,sig,u,gkx,modos,V,Lambd);
                                                  %(area,props,par,fext,glb,link,lpdva,sig,u,gkx,V,Lambd)
elseif metd==3
	[den,ugra,sigra,volgra]=sensal(area,props,par,fext,glb,link,lpdva,sig,u,vol);
    
else
    [den,ugra,sigra,volgra]=sensaj(area,props,par,fext,glb,link,lpdva,sig,u,vol);

end