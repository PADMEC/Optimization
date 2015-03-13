function [U,W_PCM]=Generate_PC_points(Mu, Cov, ft, N1)

%Comput All Permutations
    nrv = length(ft);
    if N1==1
        Permuts=ones(1,nrv);
    else
        Permuts=gridsamp([1; N1]*ones(1,nrv),N1);
    end
    npermts=N1^nrv;
    W_PCM=zeros(npermts,nrv);
    U=zeros(npermts,nrv);
    for i=1:nrv
        m=Mu(i);
        s=Cov(i,i)^.5;
        switch ft(i)
        case 1
            distype='norm';
            
        case 2
            distype='lognorm';
        case -2
            distype='lognorm';
            s(i)=-s(i);
        end 
        [xs,fs]=PC_numerc_integration(distype,m,s,N1);
        
        U(:,i)=xs(Permuts(:,i));
        W_PCM(:,i)=fs(Permuts(:,i));
    end
    
    W_PCM=prod(W_PCM,2);
end