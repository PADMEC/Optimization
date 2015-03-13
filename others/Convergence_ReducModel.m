
    irb=[ioffT:5:noffT];
    iap=0;
    for in = irb % BASIS CONVERGENCE
        iap=iap+1;
        
        D = Dt(:,1:in);
        if RBtype==1 %POD
            lamb_tol=1e-5;
             P0=mean(D,2); Ps = D - P0*ones(1,size(D,2));
             [POD_Z,nPZ,V]=PODfunction(Ps,lamb_tol);
        else %RBM
            POD_Z = D;
        end
        
        T11=cputime;
        if HYPLAS_flag
            write_Zfile(D,'RB_file.dat')
            [ rFs, rUs, rSigs, rien,rivol  ] = Hyplas_analysis( arq.x0, hyplas_command,-1, RandVarDist(2,:));
            write_Zfile(POD_Z,'RB_file.dat')
            [ pFs, pUs, pSigs, pien,pivol  ] = Hyplas_analysis( arq.x0, hyplas_command,-1, RandVarDist(2,:));
            T2=cputime;
%             [mp,ip]=max(abs((pUs-feUs)));disp(norm(pUs-feUs)/norm(feUs))
%             disp(mp)
        else
            %area=area/4;
            [pFs,pUs,pSigs,pien,pivol,pT]=NLfesol(area,props,fext,glb,tipo,RBtype,POD_Z);
            T2=cputime;
            %Test POD
%             [feFs,feUs,feSigs,feien,feivol,feT]=NLfesol(area,props,fext,glb,tipo,0);
            [rFs,rUs,rSigs,rien,rivol,rT]=NLfesol(area,props,fext,glb,tipo,RBtype,D);
            [mr,ir]=max(abs((rUs-Us)));disp(norm(rUs-Us)/norm(Us))
            [mp,ip]=max(abs((pUs-Us)));disp(norm(pUs-Us)/norm(Us))
            disp([mr,mp])
        end
        Trbm=T2-T11
%         
%         figure, plot(rUs-feUs,'.-'), hold on, plot(pUs-feUs,'o-k'), plot(feUs,'xr')

        if ioffT<noffT
            errP(iap)=abs(ien-pien);
            errR(iap)=abs(ien-rien);
            
            ap_vect=[pUs(:,end),rUs(:,end)];
            %Abs rel error
            errs=vect_err_calc(Us(:,end),ap_vect,2);
            errPdm(iap)=errs(1);
            errRdm(iap)=errs(2);
            %Norm rel error
            errs=vect_err_calc(Us(:,end),ap_vect,1);
            errPdn(iap)=errs(1);
            errRdn(iap)=errs(2);
            
%             ap_vect=[pSigs(:,end),rSigs(:,end)];
%             %Norm rel sig error
%             errs=vect_err_calc(Sigs(:,end),ap_vect,2);
%             errPs(iap)=errs(1);
%             errRs(iap)=errs(2);
            
            nPz(iap)=size(POD_Z,2);
            nRz(iap)=size(D,2);
            [errPdm(end), errRdm(end); errPdn(end), errRdn(end); nPz(end),nRz(end)]'
        end
    end %for: convergence RBM
    if ioffT<noffT
        figure, hold on
        plot(nPz,abs(errP),'r.-')
        plot(nRz,abs(errR), '.-')
        set(gca,'YScale', 'log')
        title ('Con Error')
        xlabel('N')

        figure, hold on
        plot(nPz,errPdm,'r.-')
        plot(nRz,errRdm, '.-')
        set(gca,'YScale', 'log')
        title ('Displ Max Error')
        xlabel('N')
        
        figure, hold on
        plot(nPz,errPdn,'r.-')
        plot(nRz,errRdn, '.-')
        set(gca,'YScale', 'log')
        title ('Displ Norm Error')
        xlabel('N')
        
        for i=2:8
            lamb_tol=10^-(i);
            Ps=D;
            %             P0=mean(D,2); Ps = D - P0*ones(1,size(D,2));
            %             [POD_Z,nPZ,V]=PODfunction(Ps,lamb_tol);
            [POD_Z,nPZ,V]=PODfunction(D,lamb_tol);
            nz(i)=length(nPZ);
            if HYPLAS_flag
                write_Zfile(POD_Z,'RB_file.dat')
                [pFs, pUs, pSigs, pien,pivol  ] = Hyplas_analysis( arq.x0, hyplas_command,-1, RandVarDist(2,:));
            else
                %area=area/4;
                [pFs,pUs,pSigs,pien,pivol,pT]=NLfesol(area,props,fext,glb,tipo,RBtype,POD_Z);
            end
            %Norm rel error
            err(i)=vect_err_calc(Us(:,end),pUs(:,end),1);
        end
        figure, plot(nz(2:8),err(2:8),'.-')
        set(gca,'YScale', 'log')
        title ('Displ Norm Error')
        xlabel('N')

        figure, hold on
        plot(nz(2:8),err(2:8),'.-r')
        plot(nRz,errRdn, '.-')
        set(gca,'YScale', 'log')
        ylabel('Displ Norm Error')
        legend('POD','RBM')
        xlabel('N')
   
    end
        