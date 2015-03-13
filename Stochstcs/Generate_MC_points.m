function [U,V]=Generate_MC_points(Mu, Cov, ft, n,Trans_type)

global savedState 
defaultStream = RandStream.getDefaultStream;
defaultStream.State = savedState;

        if nargin<5
            Trans_type=2;
        end
    
        Su=diag(Cov)'.^.5;

        %U = lhsnorm(Mu, Cov, n);
        Nv = length(ft);
        for i=1:Nv
            switch ft(i)
            case 1
                Unc(:,i) = lhsnorm(Mu(i), Cov(i,i), n);
            case 2
                m=Mu(i);
                s=Cov(i,i);
                MU = log(m^2 ./ sqrt(s+m^2));
                SIG = sqrt(log(s/m.^2 + 1));
                Unc(:,i) = lognrnd(MU, SIG, n,1);
            end
        end
        [par,Jvu,Juv] = preproc_randomvar(Mu,Cov, ft,Nv, Trans_type);
        V=(Unc-ones(n,1)*Mu)./(ones(n,1)*Su);
        Z=V*Jvu;
        U = Z.*(ones(n,1)*Su)+ones(n,1)*Mu;
        V = transNf(U',Juv,-ft,Nv,par);
        
end