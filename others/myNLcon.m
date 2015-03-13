 function[C,Ceq,DC,DCeq] = myNLcon(x,x0l,fext,glb,link)
%
% Fornece para o otimizador as funções restrições (normalizadas)e gradientes 
%
        global  lpdva  para   props
		global  con  i_contin xol
%         global  aprox
%         global  rest
%         global  rest2 
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
[iequa,i_contin] = compara(xol,x,itera,i_contin);

    %iequa=0;
    if iequa == 0
        [fob,gob,fre,gre] = optNLsol(props,fext,glb,link,tobj,tpres,x,ndvab,lpdva,perturb);
        xol=x;
        save tape1 x fob gob fre gre;
    else
        load tape1;
    end
% end
%
%   Pega os valores das derivadas das funcões  restrição
%
%   As restrições são normalisadas e são do tipo <=0!!
%
[C,Ceq,DC,DCeq]=trestric(fre,gre,cub,clb);
%C=max(C);
%
%	atualise para(5)  = itera 
%
	 para(5) = itera ;
     save restric C itera
