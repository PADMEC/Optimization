 function[C,Ceq,DC,DCeq] = myconrm2(x,xol,varargin)
%
% Give to the optimizer the constrainsts and their derivatives
% the functions are normalized
%
%
		 global  para 
		  global  con itp
%
%
%	pick the limits for constraints and some parameters from "para"
%
  		perturb = para(1);
		ndvab = para(2);
		tobj = para(3);
		tpres = para(4);
		itera = para(5);
%
      it = itera;

if size(con,1)>1
       clb = con(1,:);
	   cub = con(2,:);
else
    clb=[];
    cub=[];
end      
%
% Compare the design variable vectors.
% 
if find(itp==1) & sum(itp)==1 & tobj~=0
   global  dvol    
    fre=sum( x.*dvol');
    gre=dvol;
else    
  %iequa = comparabr(xol,x,itera,ndvab);
  iequa=0;
  if iequa == 0
        [fob,gob,fre,gre]=optsolrm2(tobj,tpres,x,ndvab,varargin{:});
        %t = cputime;
        %[fob,gob,fre,gre]=optsolrmEn(tobj,tpres,x,ndvab,varargin{:});
        %er = cputime-t
    xol = x ;
    %save tape1 fob gob fre gre;
  else
   load tape1; 
  end
end
%
%     pick the constraint function values 
%
%     They are normalized and <=0  type!!
%     
[C,Ceq,DC,DCeq]=trestric(fre,gre,cub,clb);
%
%	 Update para(4)  = itera 
%
	  para(5) = itera; 