function U = transNf(V,Jvu,ft,Nv,parmeters)

    %Z = Juv'*(U-Mu)';
    
    if (mean(ft)<0)
        %ft=-ft;
        Z = V';
    else
        Z = V'*Jvu;
    end
    U=zeros(size(V'));
    for i=1:Nv
        m=parmeters(i,1);
        s=parmeters(i,2);
        switch ft(i)
        case 1  % Normal distribution
            U(:,i) = Z(:,i)*s+m;
            %U(:,i) = Z(:,i)+m;
        case 2  % Lognormal distribution
            %m=Mu(i);
            %v=Cu(i,i);
            %m = log(m^2 ./ sqrt(v+m^2));
            %s = sqrt(log(v/m.^2 + 1));
            if (Z(:,i)*s+m>30)
                U(:,i) = exp(30);
            else
                U(:,i) = exp(Z(:,i)*s+m);
            end
        case -1  % inv Normal distribution
            U(:,i) = (Z(:,i)-m)/s;
        case -2  % inv Lognormal distribution
            U(:,i) = (log(Z(:,i))-m)/s;
        end
    end

    if (mean(ft)<0)
        U = U*Jvu;
    end
end
