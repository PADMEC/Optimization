
%clc,
clear

%RBM ANALYSIS
arq=arqfileGeneret();

%load  C:\MATLAB7\work\arq0


%arq.modos=TipoAnalEst;

if size(arq.x0,2)~=max(arq.arex)
    x0=arq.x0*ones(1,max(arq.arex));
else
    x0=arq.x0;
end

%rbmConverg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rbmN_Converg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SN=amostras(arq);
%hobj=gcf;
%SN=SNm(1:9,:);
%load SNrevisa
%SN=SNrev
arq=mbranalestrut(arq,SN,0,0);
% 
% if exist('plt')
%     save  C:\MATLAB7\work\arq0 arq plt
% else
%     save  C:\MATLAB7\work\arq0 arq
% end
