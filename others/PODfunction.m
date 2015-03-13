function [POD_Z,Lam,Vs]=PODfunction(D,lamb_tol)
% Proper Ortogonal Decomposition Function
%
% D - Full Snapshots Space (tryals)
% lamb_tol - principal components selection tolerance
%
% POD_Z - Principal Space Components


%Vectors Normalization
%Dmd=mean(D);
%Dnor=D - Dmd(ones(1,nDoF),:);
%Dmg=sqrt(sum(D.^2));
%Dnor=D./Dmg(ones(1,nDoF),:);

Dnor=D;

%SVD    Singular value decomposition.
[U,S,V] = svd(Dnor,0);
%D = U*S*V'. U*S=D*V. V*S'=D'*U.
%D'*U*S=D'D*V = V*S'*S=D'D*V
%V*S²=C*V, 
sig=(diag(S));
sig(sig<0)=0;
nsig=sig/sum(sig);
cumsig=cumsum(nsig);

%first k Principal Singular value selection
k=find((1-cumsig)<(lamb_tol),1);
%Dlow=U(:,1:k)*S(1:k,1:k)*V(:,1:k)';

%POD Reduced Space (reduced Order Basis)
POD_Z=U(:,1:k);%POD_Z=POD_Z(:,1:k);
Lam=nsig(1:k);

Sm1=diag(1./diag(S));% Sm1 = S^-1
Vs = V*Sm1(:,1:k);
%figure, plot(pv-U(:,1:k),'.-')

%% Eig-Test
% %Covariance Matrix (test)
% Cd=Dnor'*Dnor;
%  
% % %Eigenvectors decomposition
% [Vs,Lam]=eig(Cd);
% norm((Cd*Vs)-(Vs*Lam))
% norm((Cd*V)-(V*S.^2))
% tolp=
% [nDoF]=size(D);% 
% figure, plot(real(log(diag(S.^2))),'.-')
% hold on,plot(real(log(diag(Lam(nDoF(2):-1:1,nDoF(2):-1:1)))),'o--r')
% 
% Lam=diag(Lam)';
% Lam(Lam<0)=0;
% Lam=realsqrt(Lam(k:-1:1));
% Vs=Vs(:,(k:-1:1));
% nLam=Lam/sum(Lam);

% %Principal Components Selection
% totL=1;ilam=1;
% while totL>lamb_tol
%     totL=totL-nLam(ilam);
%     ilam=ilam+1;
% end
% ilam=ilam-1;
%Lam(ilam:nZ)

%POD Space
%POD_Z=Dnor*Vs(:,ilam:nZ);
%POD_Z=Dnor'*Vs(:,ilam:nDoF);

% POD_Z2=Vs(:,ilam:nDoF);
% POD_Z=Ul;

%POD_Z2=Dnor*(Vs(:,1:ilam)./Lam(ones(nZ,1),1:ilam));

%APROXIMACAO!!!
%Dn=~POD*a + Dm;

%%TEST
% a=POD_Z\Dnor(:,1);
% a2=POD_Z2\Dnor(:,1);


%POD_Z=Dlow;