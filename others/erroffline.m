function[co,tjjqq,vjq]  =  erroffline(Z)
%
%
%  offline stage for error analysis:
%
%
   load trussdata.mat
% computes khat matrix
%
N=size(Z,2);
Khat=K1{1}-K1{1};
for i=1:size(K1,2)
    Khat = Khat+K1{i};
end
%   calculates vectors zq 

%
%  solves for each dv in turn
%
     jcol = 0;
     for idvab = 1:size(K1,2)
%       
          for in = 1:N
%
%              get solutions for each of the sampling in 
%
               x  = Z(:,in);
               fq = K1{idvab}*x;
               zq = full(Khat\fq);
%
%              builds the matrix zhat
%
               jcol = jcol + 1;
               zhat(:,jcol) = zq;
           end
      end
      zh = zhat;
%
%     computes zohat, co
%
      zohat =  full(Khat\F);
      co = zohat'*F;
%
%      computes vjq vector
%
      vjq = zhat'*F;
%
%      computes Tjj'qq' matrix
%      
      for ii = 1:size(K1,2)*N
          for jj = 1:size(K1,2)*N
              y = zhat(:,ii);
              x = zhat(:,jj);
              tjjqq(jj,ii) = y'*Khat*x;
          end
      end
    
     
 
  