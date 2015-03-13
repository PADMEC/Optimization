function [C,Ceq,DC,DCeq]=trestric(fre,gre,cub,clb)


  %nsize = size(cub,1);
  dv=abs(cub)+abs(clb);dv(dv==0)=1;
  redv=isinf(cub)|isinf(clb);
  clb(isinf(clb))=1e12;
  cub(isinf(cub))=1e12;
  dv(redv)=1;%cub(cub==inf)=1e12;
  if length(fre)==length(cub)
      rnor=(fre'-cub)./dv;
      rnor=[rnor (clb-fre')./dv];
      
  elseif length(fre)>length(cub)
      lc=length(cub);
      if lc==0, rnor=fre';
      else
          dl=length(fre)-lc;
          rnor=(fre(1:lc)'-cub)./dv;
          rnor=[rnor (clb-fre(1:lc)')./dv];
          rnor=[rnor fre(lc+1:lc+dl)'];
      end
  end
  rnor(isnan(rnor))=-1;
  rnor(rnor==-inf)=-1e12;
   % nsize = size(fre,1) ;
   % for isize = 1:nsize 
   %     icon = isize ;
   %     rnor(icon)         = (fre(isize) - cub(isize))/cub(isize);
   %     rnor(icon + nsize) = (clb(isize) - fre(isize))/abs(clb(isize));
   % end
%
%   finalmente as funcões restrição 

    C = rnor ;
    Ceq=[];
%
%
%   compute das derivadas das restrições  normalisadas:
%
[ndvab,nsize] = size(gre)	;
gnor=[];
if length(cub)==nsize
    for i=1:ndvab
      gnub(i,:)=gre(i,:)./dv;
      gnlb(i,:)=-gre(i,:)./dv;
    end
      gnor = [gnub gnlb];
  elseif length(cub)<nsize
      lc=length(cub);
      if lc==0, gnor=gre;
      else
          for i=1:ndvab
              gnub(i,:)=gre(i,1:lc)./dv;
              gnlb(i,:)=-gre(i,1:lc)./dv;
          end
          gnor = [gnub gnlb gre(:,lc+1:end)];
      end
end
      
%   1. inicialise:
%
    %DC = zeros(ndvab,2*nsize);
%
%   2. atualise:
%
    DC = gnor;
    DCeq=[];