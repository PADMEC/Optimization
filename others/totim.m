function [fre,gre,fob,gob]=totim(tres,tobj,vol,dvol,sn,dsn,sig,dsig,u,du,floc,gfloc,ndvab,x,props,lamb,dlamb)

    global  id mo fs mf nor beta l itp

        global Messac
        
    fre=[];frei=[];
    gre=[];grei=[];

    for i=itp
        if i==1             %VOLUME
            frei = vol;
            grei = dvol;
            
        elseif i==2         %FLAMBAGEM GLOBAL
            frei = lamb;
            grei = dlamb;

        elseif i==3         %FLAMBAGEM LOCAL
            frei =floc;
            grei=gfloc;
            
        elseif i==4         %TENSAO
            frei = sig;
            grei = dsig;
            
        elseif i==5         %DESLOCAMENTO
            frei = u;
            grei = du;

        elseif i == 6        %TENS  MAX
            
%             global SIGu
%             
%             [ms,imd]=max(abs(sig));
%             %gms = dsig(:,imd)*sign(sig(imd))/ms;
%             gms = dsig(:,imd)*sign(sig(imd));
%             
% 
%             sig(imd)=[];
%             dsig(:,imd)=[];
%             %fra=abs(sig)/ms-1;
%             fra=(abs(sig)-ms)/SIGu;
%             for is=1:length(sig)
%                 %gra(:,is)= dsig(:,is)*sign(sig(is))/ms;
%                 gra(:,is)= (dsig(:,is)*sign(sig(is)) - gms)/SIGu;
%                 % % d(s/ms)/dx = gs/ms - s/(ms^2)*gms = (ds*ms - s*dms)/(ms^2) 
%                 %ds = dsig(:,is)*sign(sig(is))*ms;
%                 %dms = sig(is)*gms;
%                 %gra(:,is)= (ds - dms)/(ms^2);
% 
%             end
%             fre=[fre; fra];
%             gre=[gre gra];
            
        elseif i==7         %DESLOCAMENTO ESPECIFICO
            frei = u(id);
            grei = du(:,id);
        elseif i==8         %TENSAO ESPECIFICA
            frei = sig(id);
            grei = dsig(:,id);

        end
        fre=[fre;frei];
        gre=[gre grei];
    end
    
    if tobj == 1            %VOLUME
        fob = vol;
        gob = dvol;
    elseif tobj == 2        %ENERGia
        fob = sn;
        gob = dsn;
    elseif tobj==3          %FLAMBAGEM GLOBAL
            fob = -lamb(1);
            gob = -dlamb(:,1);
        
    elseif tobj == 4        %NORMA DOS DESLCAMENTO
        fob = norm(u);
        C=sum(u.^2)^(-.5)/2;
        for idvab = 1:ndvab
            gob(idvab) = C*2*du(idvab,:)*u;
        end
        
    elseif tobj == 5        %DESL.  MAXIMO
        [fob,imd]=max(abs(u));
        gob = du(:,imd)*sign(u(imd));
        %rdu=du(:,idm)*sign(u(idm));
        u(imd)=[];
        du(:,imd)=[];
        fra=abs(u)-fob;
        for i=1:length(u)
            gra(:,i)=du(:,i)*sign(u(i))-gob;
        end
        fre=[fre; fra];
        gre=[gre gra];


    elseif tobj == 6        %TENS  MAX
        [fob,imd]=max(abs(sig));
        gob = dsig(:,imd)*sign(sig(imd));
        
        sig(imd)=[];
        dsig(:,imd)=[];
        fra=abs(sig)-fob;
        for i=1:length(sig)
            gra(:,i)=dsig(:,i)*sign(sig(i))-gob;
        end
        fre=[fre; fra];
        gre=[gre gra];
        %imd=3;
        %fob=abs(sig(imd));
        %gob = dsig(:,imd)*sign(sig(imd)));

    elseif tobj == 7        %DESL.  ESP
        
        global kb
        global idm
        
        %fob = [u(idm)^2/2' sig(kb)^2/2'];
        fob = abs([u(idm)' sig(kb)']);
        if isempty(kb)
            %rdu=du(:,idm)*u(idm);
            rdu=du(:,idm)*sign(u(idm));
            rdt=[];
        elseif isempty(idm)
            rdu=[];
            %rdt=dsig(:,kb)*sig(kb)';
            rdt=dsig(:,kb)*sign(sig(kb))';
        end

        gob = sum([rdu rdt],2);
        
        if ~isempty(Messac)
            %%%%%%%% Exemplo A. Messac %%%%%%
        
            fob = sum(abs(u).*[.75 .25]');
            rdu=[];
                for i = 1:ndvab
                    rdu(i,:) = du(i,1:2).*sign(u(1:2))'.*[.75 .25];
                end
            gob = sum([rdu],2);
        end
        %%%%%%%% Exemplo A. Messac %%%%%%
            %%%%%%%%  FIM  %%%%%%
        
    elseif tobj == 0%MULTI-OBJETIVO
        %%%%%%%%    MULTI-OBJETIVO  %%%%%%%%
        
            %ENERGIA - 2
        foben = sn;
        den=dsn;
            
            %norma DESL. - 4
        if mo(1,4)==1
            fobdt = norm(u);
            C=sum(u.^2)^(-.5)/2;
            for idvab = 1:ndvab
                rdu(idvab,:) = C*2*du(idvab,:)*u;
            end
            ddt = rdu;
            else fobdt=0;ddt=dvol-dvol;
        end
            
            %DESL. MAX - 5
        [fobdm,imd]=max(abs(u));
        ddm = du(:,imd).*abs(u(imd))./u(imd);

            %TENS. MAX - 6
            %imd=3;
        [fobtm,imd]=max(abs(sig));
        dtm = dsig(:,imd).*sign(sig(imd));

            %DESL. ou Tensao ESP. - 7
        if mo(1,end)>0
            global kb
            global idm
            %fobde = [u(idm)'.^2/2 sig(kb)'.^2/2];
            fobde = [abs(u(idm)') abs(sig(kb)')];
            if isempty(kb)
                for i = 1:ndvab
                    %rdu(i,:) = du(i,idm).*u(idm)';
                    rdu(i,:) = du(i,idm).*sign(u(idm)');
                end
                rdt=[];
            elseif isempty(idm)
                rdu=[];
                for i = 1:ndvab
                    %rdt(i,:) = dsig(i,kb).*sig(kb)';
                    rdt(i,:) = dsig(i,kb).*sign(sig(kb))';
                end
            end
            dde = [rdu rdt];
            
        if ~isempty(Messac)
                %%%%%% Exemplo A. Messac %%%%%%
                    dde=[];
                    fobde =[];
                    rdu=[];
                    fobde = sum(abs(u).*[.75 .25]');
                        for i = 1:ndvab
                            rdu(i,:) = du(i,1:2).*sign(u(1:2))'.*[.75 .25];
                        end
                    dde = sum([rdu],2);
                %%%%%% Exemplo A. Messac %%%%%%
        end
            %%%%%%  FIM  %%%%%%        
            
        else
            fobde=0;
            dde=dvol-dvol;
        end
        
        f=[vol foben -lamb(1) fobdt fobdm fobtm fobde];
        df=[dvol den' -dlamb(:,1) ddt ddm dtm dde];
        imo=find(mo(1,:));
        
        
%%%%%%%%%%%%%%%% FUNÇAO SUBSTITUTA
        if fs==1    %f=f(x)/f0, F=sum(f);
            %Soma Ponderada
            
            minf=abs(min(mo(3:end,imo)));
            
            fob=sum(f(imo).*beta./minf);
            
            ddfx=df(:,imo);
            for i=1:ndvab
                gob(i)=sum(ddfx(i,:).*beta./minf);
            end
            
        elseif fs==2    %f=(f(x)-minf)/(maxf-minf), F=max(f)
            % Min-Max ??????????????????????????            
            
            %% NNC %%
            fobs=(f(imo)-mf)./l.*beta;
            df=df(:,imo);
            for i=1:length(imo)
                gobs(:,i)=df(:,i)*beta(i)/l(i);
            end
            [fob,ifo]=max(fobs);
            gob=df(:,ifo)*beta(ifo)/l(ifo);
            %[fob,ifo]=max(dfx.*beta./dfo);
            %ddfx=df(:,imo);
            %gob=ddfx(:,ifo)'*beta(ifo)/df(ifo);
            
            fobs(ifo)=[];
            gobs(:,ifo)=[];
            fra=fobs-fob;
            for i=1:length(fobs)
                gra(:,i)=gobs(:,i)-gob;
            end
            fre=[fre; fra'];
            gre=[gre gra];

            
        elseif fs==3|fs==5 % NBI
            
            
            f=f(imo)';
            df=df(:,imo);
            
            if fs==31
                %F=[1;f-mf];
                b=Am1*F;
                fob=-b(end);
                for i=1:ndvab
                    gf(i)=-Am1(end,2:end)*(df(:,i));
                end
                gob=gf;
            end
            global t

            fob=t;
            gob=[zeros(1,ndvab) 1];
            %res=max(t)-min(t);
            if fs==5
                global par
                res=beta-t*nor(par)-f(par)./mf';
                for i=1:ndvab
                    gres(i,:)=-df(i,par)./mf;
                end
                gre = [gre           gres;
                zeros(1,size(gre,2)) -nor(par)'];
            else
                res=beta-t*nor-f./mf;
                for i=1:ndvab
                    gres(i,:)=-df(i,:)./mf';
                end
                gre = [gre           gres;
                zeros(1,size(gre,2)) -nor'];
            end
            fre = [fre; res];
            
        elseif fs==6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            f=f(imo)';
            fob=f(1);
            df=df(:,imo);
            gob=df(:,1);
            res=f(2:end)./beta(2:end)'-1;
            for i=1:ndvab
                gres(i,:)=df(i,2:end)./beta(2:end);
            end
            fre = [fre; res];
            gre = [gre  gres];

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        elseif fs==4 % NNC
            
            %fob=(f(max(imo))-mf(end))./l(end);
            %gob=df(:,max(imo))./l(end);
            fn=(f(imo)-mf)./l;
            dfn=df(:,imo);
            for i=1:ndvab
               dfn(i,:)=dfn(i,:)./l;
            end
            fob=fn(end);
            gob=dfn(:,end);
            
            gres=dfn*nor';
            res=nor*(fn-beta)';
            fre = [fre; res];
            gre = [gre  gres];
            
        end
        
    end

    