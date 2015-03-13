function [iequa,i_contin] = compara(olx,x,itera,i_contin)
%
%	Compara valores das variáveis de projeto p/ identificar se é
%	 				de fato um novo projeto.
%
%
%
%
if itera ==0
    iequa = 0;
    i_contin = 0;
else
    epsaf = (eps)*200;
    iequa = 1 ;
    if i_contin ~= 2, i_contin = 1;end
    lx=abs(x); lx(x==0)=1;
    xstot = abs(olx - x) ./ lx;
    if any(xstot>= epsaf)
        iequa = 0;
    end
    if any(xstot>= 1e-2)
        i_contin = 0;
    end
end
