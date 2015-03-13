
function [parmeters, J, J1,Sig] = preproc_randomvar(Mu,Cu, ft,Nv, Trans_type)


    Su=diag(Cu).^.5;
    Sig=diag(1./Su);
    C=Cu./(Su*Su');
    [J,J1] = JacobCorMat(C,Trans_type);
    
    % Parameters of the pdfs - probabilite distribuction fuction
    parmeters=zeros(Nv,5);
    for i=1:Nv
        switch ft(i)
        case 1  % Normal distribution
            m=Mu(i);s=Cu(i,i)^.5;
            parmeters(i,1:2)=[m,s];
            
        case 2  % Lognormal distribution
            m=Mu(i);v=Cu(i,i);
            mu = log(m^2 ./ sqrt(v+m^2));
            s = sqrt(log(v/(m.^2) + 1));
            
            parmeters(i,1:2)=[mu,s];
            
        end
    end
    
    
end


function [A,A1] = JacobCorMat(Cu,Trans_type)

% V = A'U, v eh nao-correlacionada
% Cv = cov(V,V') = cov(A'U,U'A) = A'cov(U,U')A
% Cv = A'Cu.A -> Diag
% A.Cv = Cu.A
% [A, D] = eig(Cu)
% A1 = A*(D^-.5) -> covariancia unitaria
% Cholesky
% A1'Cu.A1 = I -> Cu = L'.L
% A1 = L
% Cu = (A1'^-1) * (A1^-1)
% Cu = L'*L (Cholesky)

switch Trans_type
case 1
    A = chol(Cu);
case 2
    [A, D] = eig(Cu);
    A = (A*(D^-.5))^-1;

case 3
    [A, D] = eig(Cu);
    A=A';
end
A1=A^-1;

% V = lhsnorm(zeros(2,1), eye(2), 1000);
% U = V*A';
% hold on, plot(U(:,1),U(:,2),'.g');
end