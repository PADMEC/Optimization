 function[F,G] = myfrm2(x,xol,varargin) 
% Give to the optimizer the constrainsts and their derivatives
% the functions are normalized
%
%
%        global  kns
		 global  para 
		  global itp
          global icount
          global xt
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
     % clb = con(1,:);
	 %  cub = con(2,:);
      it = itera;
     
%
% Compare the design variable vectors.
% 
if tobj==1
    global  dvol
    fob=sum( x.*dvol');
    gob=dvol;
else
   %iequa = comparabr(xol,x,itera,ndvab);
   itera = itera + 1 ;
   iequa=0;
  if iequa == 0
    icount = icount + 1;
    [fob,gob,fre,gre]=optsolrm2(tobj,tpres,x,ndvab,varargin{:});
    %t = cputime;
    %[fob,gob,fre,gre]=optsolrmEn(tobj,tpres,x,ndvab,varargin{:});
    %ef = cputime-t
    xol = x ;
    if itera > 1
           %xt(itera,:) = xol;
        itt = itera;
    end     
    %save tape1 fob gob fre gre;     
  else 
   load tape1;
  end
end
%
%   Picks the objective function value and its derivatives
%
%    objective function:
%
    F = fob;
%
%   Gradient:

    G = gob;

%	 atualize para(4)  = itera 
%
	  para(5) = itera ;
	  itfun = itera	;
