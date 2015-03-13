function Jzv = inv_normeq(Mu,Su,U,ft,pdfun)

    nv=length(Mu);
    for i=1:nv
        switch ft(i)
        case 1  % Normal distribution
            Jzv(i,i) = Su(i);
        case 2  % Lognormal distribution
            ksi=sqrt( log( 1 + ( Su(i) ./ Mu(i) ).^2 ) ).*U(i);
            Jzv(i,i)=ksi;
        case 3  % other distribution
           %lambda = marg(i,5);
           %k = marg(i,6);
           %pdf1 = lambda * (lambda*x(i))^(k-1) / gamma(k) * exp(-lambda*x(i));
           pdf1 = feval(pdfun,U(i),Mu(i),Su(i));
           z=(U(i)-Mu(i))/Su(i);
           pdf2 = normpdf(z);
           Jzv(i,i) = pdf2/pdf1;
        end
    end
    
end