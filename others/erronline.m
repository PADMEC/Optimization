function[dn]  =  erronline(alph,co,munew,ndvab,tjjqq,vjq)
%
%
%  online stage for error analysis:
%
%
vjq;
%
%   picks the smallest value of munew
%
  g = min(munew);
%   
% initialize
%  
   pt2 = 0;
   pt3 = 0;
%   calculates part 2 of error equation first!!
% atencao (particular case ndvab = 3 )
%
       
%
%  solves for each dv in turn
%
N=length(alph);
     jcol = 0;
     for idvab = 1:ndvab;     
          for in = 1:N;
              jcol = jcol + 1;
              cvjq = vjq(jcol);
              pt2 = pt2 + alph(in)*munew(idvab)*cvjq;
%
%              Then calculates part 3 of error equation 
%
               irow = 0;
               for jdvab = 1:ndvab; 
                   for jn = 1:N;
                       irow = irow + 1;
                       ctjjqq = tjjqq(irow,jcol);
                       pt3 = pt3 + alpha(in)*alph(jn)*munew(idvab)*munew(jdvab)*ctjjqq;
                   end
               end
           end
       end
             co;
             pt2;
             pt3;

%
%     finally computes the error dn:
%
     dn = (1/g)*(co +2*pt2 + pt3);
    