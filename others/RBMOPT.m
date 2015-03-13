


clc,clear

% RUN RBM OPTIMIZATION
%arq=arqfileGeneret();

load arq0

%arq=tipotim(arq);

%SN=amostra(arq)

ndv=max(arq.arex);
xmin=arq.xmin*ones(1,ndv);
xmax=arq.xmax*ones(1,ndv);
SN=amostraRB(ndv,xmin,xmax,6);

%load amostra_12_leo.mat
%SN=S12;
%disp('BASES :')
%SN
% hobj=gcf;

plt='';
arq=mbr(arq,plt,SN);
%save  C:\MATLAB7\work\arq0 arq plt

