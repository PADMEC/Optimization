 function[F,G] = myNLfun(x,x0l,fext,glb,link) 
%
% Fornece para o otimizador a função objetivo e gradiente 
%
%   globals
%   
        global  lpdva  para  props
        global  i_contin xol Xs Fobs
%           global  apr
%           global  objet
%
%	pega  nas variaveis  "para"  , "props"  e "con" alguns paramentros adicionais
%
    perturb = para(1);
	ndvab   = para(2);
	tobj    = para(3);
	tpres   = para(4);
	itera   = para(5);

%
%
% 
% Compara os vetores das variáveis de projeto
% 
% if tobj==1
%     comp=props(2,:);
%     ro=props(4,:);
%     fob=sum( x(link(:,1)') .* comp.*ro );
%     for i=1:ndvab
%         bardv=find(link(:,1)==i);
%         gob(i,1)=sum(comp(bardv).*ro(bardv));
%     end
% else
    [iequa,i_contin] = compara(xol,x,itera,i_contin);
    itera = itera + 1 ;
    %iequa=0;
%     xi = xol;
%     xf = x;

    if iequa == 0
        [fob,gob,fre,gre]=optNLsol(props,fext,glb,link,tobj,tpres,x,ndvab,lpdva,perturb);
	    xol = x ;
        save tape1 fob gob fre gre x;
        
    else
        load tape1;
    end
% end
   % Xs=[Xs; xol];
   % Fobs=[Fobs; fob];
%
%   Pega os valores das funcões objetivo e sua derivada 
%
%   funcão objetivo:
%
    F = fob;
%
%   Gradiente:
%
    G = gob;
    
% func = [ F G ]
%
%   atualise para(5)  = itera 
%
	para(5) = itera ;
% 	itfun   = itera ;
    