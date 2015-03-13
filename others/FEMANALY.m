

%clc,
%clear('all')
clear variables global

%profile on, profile clear

%FEM ANALYSIS
%arq=arqfileGeneret();
%load ExempA-Messac.mat
%load 200barLeoRB.mat
%load resultsOtimo.mat
% 
% load arq960bar.mat
% arq.x0=0.05;
% 
% load arq64bar.mat
% arq.x0=0.05;
% 
% load 1800EXE(LargTruSistm).mat
% k_p =1/13;

 load  NewPortcLarg
 arq.x0=1e-3;
 arq.E=2.07e7;
 k_p =10000; % Relationship between FxMatCurve
 nOff=24; %n of non-lin analyses offline
 Lx=2; %up bound of x
 
arq.x0=[.1 .6 .8];
% - New Portc opt design variables
Xopt=[ 0.2172    1.0279    1.3598];
arq.x0=Xopt;
 
 plotlogc=1;
 
% arq.x0=0.05;
% 
% load MalhaOtimPonte.mat
% arq.x0=arq.x;

%load NonLin_T.mat

%cm, N 1KN/m2=1000/10^4(N/cm2)=10(N/cm2)
%load Flamb3X.mat

arq=MatProp(arq);


[arq,area,props,fext,glb,tipo,nel,mdf,mdv]=pre_compus_arq(arq);
 fext =fext*k_p;


%%%%%%%%%%%%%%%%%%% Failure Analises ?%%%%%%%%%%%%%%%%%%%%%

global es0
es0=zeros(nel,1);
%%%Simple FEM Analysis
%[Fs,Us,Sigs,arq.ien,arq.ivol]=NLfesol(area,props,fext,glb,tipo);

[D,Fx,Zkn,Cs,Toff,x0s,k_iter]=OfflinePOD(arq,props,fext,glb,tipo,nOff,Lx);

%Kd=f
%f0+df/du.du=f, f-fo=R
%
%dF.du=R
%dF.(Z.a)=R
%Z'.dF.Z=Z'.R
nZ=size(D,2);
area = arq.x0(arq.arex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ELASTIC LINEAR TEST %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[lu,lsig,len,vol,floc]=fesolOPT(area,props,fext,glb,tipo);
%[rblu,rblsig,rblen,vol,floc]=fesolOPT(area,props,fext,glb,tipo,1,D);

%%%Simple FEM Analysis
T11=cputime;
disp('FEM Analyses')
[Fs,Us,Sigs,arq.ien,arq.ivol,T]=NLfesol(area,props,fext,glb,tipo);
T2=cputime;
Tfem=T2-T11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% RBM - Reduced Basis Method %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T21=cputime;
% [Fs,Us,Sigs,arq.ien,arq.ivol]=NLfesol(area,props,fext,glb,tipo,1,D,Zdu);
% T3=cputime;
% Trbm=T3-T21
% 
% plot_NL_Results(Fs,Us,Sigs,'RBM',id_plot)
% plot_sig_vet(Fs(:,end),0*Us(:,end),Sigs(:,end),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% POD - Proper Ortogonal Decomposition %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T31=cputime;

lamb_tol=1e-5;

%POD_Z=eye(1100);

POD_Z=D;
[POD_Z,nPZ,V]=PODfunction(D,lamb_tol);
% % POD_Z = D*V
% 
% nitt=size(D,2);
% for ipz=1:nPZ
%     for iit=1:nitt
%         ZK{ipz}=ZK{ipz}+Zkn{iit}*V(ipz,iit);
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPACE TEST 
% Z=D;
% Z=POD_Z;
% ai=Z\Us;
% di=Z*ai-Us;
% ndif=(sum(di.^2).^.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T32=cputime;
disp('POD Analyses')
%[Fs2,Us2,Sigs2,arq.ien2,arq.ivol2]=Lfesol(area,props,fext,glb,tipo,1,POD_Z,POD_Zdu);
[Fs2,Us2,Sigs2,arq.ien2,arq.ivol2,T2]=NLfesol(area,props,fext,glb,tipo,1,POD_Z,Zkn,Cs,arq.x0);
%[Fs,Us,Sigs,arq.ien,arq.ivol]=NLfesol(area,props,fext,glb,tipo);
T4=cputime;

Tpre_pod=T32-T31
Tpod=T4-T32

%profile report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotlogc
ms=max(abs(Sigs));ms2=max(abs(Sigs2));
md=max(abs(Us));md2=max(abs(Us2));
mf=max(abs(Fs));mf2=max(abs(Fs2));
figure, plot(md,ms,'.-'),hold on, plot(md2,ms2,'.-r')
xlabel('max u'), ylabel('max \sigma'),legend('FEM','POD')
title(['POD - (N e tol): ' num2str([size(POD_Z,2) lamb_tol])])

%figure, plot(T-T11,mf,'.-'),hold on, plot(T2-T32,mf2,'.-r')
figure, plot(md,T-T11,'.-'),hold on, plot(md2,T2-T32,'.-r')
ylabel('Time'), xlabel('max d'),legend('FEM','POD')
title(['POD - (N e tol): ' num2str([size(POD_Z,2) lamb_tol])])

figure, plot(mf,ms,'.-'),hold on, plot(mf2,ms2,'.-r')
xlabel('F'), ylabel('max \sigma'),legend('FEM','POD')
title(['POD - (N e tol): ' num2str([size(POD_Z,2) lamb_tol])])
%id_plot=plot_NL_Results(Fs,Us,Sigs,'FEM');
% figure,plot_sig_vet(area,Fs(:,end),Us(:,end),Sigs(:,end),1);
% aviname='NonLin_FEM';type=1;
% PlasticMove(Fs,Us,Sigs,type,aviname)

%plot_NL_Results(Fs2,Us2,Sigs2,'POD',id_plot)
% figure,plot_sig_vet(area,Fs2(:,end),Us2(:,end),Sigs2(:,end),2);
%aviname='NonLin_POD';type=1;
%PlasticMove(Fs,Us,Sigs,type,aviname)
end
