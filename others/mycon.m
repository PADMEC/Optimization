 function[C,Ceq,DC,DCeq] = mycon(x,xol,fext,glb,link)
%
% Fornece para o otimizador as funções restrições (normalizadas)e gradientes 
%
        global  lpdva
		global  para 
		global  props
		global  con itp
        global  aprox
        global  rest
        global  rest2 
%
%	pega  nas variaveis  "para"   "props" e "con" alguns paramentros adicionais
%
    perturb = para(1);
	ndvab = para(2);  	   
	tobj = para(3); 	  
	tpres = para(4);  
	itera = para(5);  
%
if size(con,1)>1
       clb = con(1,:);
	   cub = con(2,:);
else
    clb=[];
    cub=[];
end    
%
% Compara os vetores das variáveis de projeto
% 
iequa = compara(xol,x,itera,ndvab);
if find(itp==1) & sum(itp)==1 & tobj~=0
    comp=props(2,:);
    ro=props(4,:);
    fre=sum( x(link(:,1)) .* comp.*ro );
    for i=1:ndvab
        bardv=find(link(:,1)==i);
        gre(i,1)=sum(comp(bardv).*ro(bardv));
    end
else    
    %iequa=0;
    if iequa == 0
        [fob,gob,fre,gre] = optsol(props,fext,glb,link,tobj,tpres,x,ndvab,lpdva,perturb);
        save tape1 x fob gob fre gre;
    else
        load tape1;
    end
end
%
%   Pega os valores das derivadas das funcões  restrição
%
%   As restrições são normalisadas e são do tipo <=0!!
%
[C,Ceq,DC,DCeq]=trestric(fre,gre,cub,clb);
%
%	atualise para(5)  = itera 
%
	 para(5) = itera ;
     save restric C itera

% novo !!!

% 1 linear; 2 reciproca; 3 rec modificada; 4 conservativa; 5 convexa; 6 polinomial
% qual o tipo de aproximacao
%iapr = 2 ;

%xi = aprox(1);
%xf = aprox(2);

%[freap,greap] = aprox (iapr,fre,gre,xi,xf,para);


%format short e;
%rest = [ C DC ]
%rest2 = [ freap greap ]
%format;