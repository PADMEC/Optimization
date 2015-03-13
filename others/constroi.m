% sub-rotina para construir a matriz LD 

%---------------------------------------------------------
% nnos = numero de nos 
% ngl = numero de graus de liberdade
% ID entra como matriz de zeros e uns e sai com uma matriz
% contendo os indices dos graus de liberdade e
% zeros nas restriçoes  
% con = contem a conectividade de cada elemento 
% --------------------------------------------------------
function LD = constroi(ID,conect)

[nnos,n] = size(ID);
ngl = 0.0;

for i = 1:nnos
	for j = 1:n
		g = ID(i,j);
		if g >= 1
			ID(i,j)=0.0;
		else
			ngl = ngl + 1;
			ID(i,j) = ngl;
		end;
	end;
end;

[nelm] = length(conect);
for i = 1:nelm,
	for j = 1:n,
		LD(i,j) = ID(conect(i,1),j);
		LD(i,j+n) = ID(conect(i,2),j);
	end;
end;

LD = LD';