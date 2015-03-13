global echoReliabPlot i_contin grad_calculation V_pre_optimum k_LS

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FORM
        %
        %U=DinamcVar;% First guess
        %DinamcVar=0;
        %DinamcVar=1;% Dual procedure
        if strcmp(Method,'R2DO')
            DinamcVar=0;%i_contin;
        end
        if DinamcVar
            % Set FORM start  point
            if isempty(V_pre_optimum)
                % no FORM start point 
                U=Mu;
                V = U'*0;
                iscont=0;
                %U=DinamcVar;% First guess
                %DinamcVar=0;% Dual procedure
            else
                % V_pre_optimum FORM start point 
                V = V_pre_optimum;
                U = transNf(V,Jzu,ft,Nv,parmeters);
                iscont=1;
            end
        else
            U=Mu;
            V = U'*0;
            iscont=0;
            %V = transNf(U,Juz,-ft,Nv,parmeters);
        end
        
        % COMPUTE G FIRST POINT
        [gV,dgU] = feval(G_fun,U,varargin{:});
        if DinamcVar
            gV=-1;
        end
        sgnG=sign(gV);

        if ((gV > 0) && (~iscont))
            % Insigth!!!!!
            % Too failure region, aproximate via normal distribuction
            beta = -gV*10;
            i=1;
            
            Us{1}(1,:)=U;
            Vs{1}(1,:)=V';
            %dUs(1,:)=U;
            %dVs(1,:)=V';
            gs(1,1)=gV;

            % dbeta = -dgx/|dgV| -> -dgx
            %dgV = dgU/norm(dgU);
            
            grad_calculation = 1;
            var_type=1;
            [Gs,dgx,varargout] = feval(G_fun,U,varargin{:});
            var_type=0;
            grad_calculation = 0;
            
            % Computing: dGv/dV 
            Jzv = inv_normeq(Mu,Su,U,ft);
            dgV = Jzu*Jzv*dgU;
            dbeta = -dgx/norm(dgV)*10;
            PF = normcdf(-beta );
            
            fprintf('Infesible, Norm G constraint: %g, beta = %g \n', max(gV),norm(V))

        else

            if echoReliabPlot
                figure,hold on, plot3(U(1),U(2), gV,'or')
            end

            % Computing: dGv/dV 
            Jzv = inv_normeq(Mu,Su,U,ft);
            dgV = Jzu*Jzv*dgU;

            normg=norm(gV);
            i=0;
            Us=U;
            Vs=V';
            gs=gV;

            if echoReliab
                fprintf('Iter %d, Norm G constraint: %g, beta = %g \n', i, gV,norm(V))
            end
            
            % begin MPP search
            while (normg>1e-4)
            %% Begin FORM Iterations
                %dgU = FiniteDiff(G_fun,U,gV,varargin{:});
                % dG/dV = dG/dU*dU/dV = dG/dU*Jvu
    %             du = -gV*dgU/(dgU'*dgU);
    %             U = U + du';

                Jzv = inv_normeq(Mu,Su,U,ft);

                %dg/dv=dg/du*du/dv, du/dv=Jvu
                dgV = Jzu*Jzv*dgU;
                V0=V;

                %% HL-RF
                if norm(dgV)==0
                    V = (dgV'*V - gV);
                end
                V = (dgV'*V - gV)*dgV/(dgV'*dgV);
                dV = V-V0;

                %% Next step
                U0=U;dgU0=dgU;gV0=gV;dgV0=dgV;
                U = transNf(V,Jzu,ft,Nv,parmeters);
                [gV,dgU] = feval(G_fun,U,varargin{:});
                Jzv = inv_normeq(Mu,Su,U,ft);
                dgV = Jzu*Jzv*dgU;

                
                if echoReliabPlot
                    plot3(U(1),U(2), gV,'.r')
                end

                %% Busca linear
                while (gV0<-0.001)&&(gV>0.001)
                    if echoReliab
                        fprintf('Busca lin: %d, Norm G constraint: %g, beta = %g \n', i, gV,norm(V))
                    end

                    i=i+1;
                    Us(i+1,:)=U;
                    Vs(i+1,:)=V';
                    gs(i+1)=gV;

                    %APPROXIMATIONS
                    %
                    %Linear: 
                    %(gV-gV0)*h+gV0=0
                    h_lin=gV0/(gV0-gV);
                    %Nonlinear (quadratic) starting from init
                    % dgU0 + Hes*dU = dgU, Hes.dU = ddgU
                    % gV0 + dgU*dU*h + dU'H*dU/2*h^2 = 0
                    % gV0 + dgU*dU*h + dU'ddgU/2*h^2 = 0
                    % gV0 + dgU*dU*h + dU'ddgU/2*h^2 = 0
                    % a = dU'ddgU/2, b = dgU*dU, c = gV0;
                    dU = U-U0;
                    ddgU = dgU - dgU0;
                    dgU2 = (dgU+dgU0)/2; gV2=(gV+gV0)/2;
                    a = dU*ddgU/2;
                    b = dU*dgU2;
                    c = gV0;
                    %x^2 + 2bx + c = (x + b)^2 = - c + b^2
                    %(x + b)^2 =  b^2 - c, x = (b^2-c)^.5 - b;
                    if abs(a)<1e-8
                        if abs(b)<1e-8
                            h=0;
                        else
                            h=-c/b;
                        end
                    else
                        rs=roots([a,b,c]);
                        h = min (rs(rs>-.5));
                        if isempty(h)
                            h=-b/2/a;
                        end
                        if ~isreal(h)
                            h=-b/2/a;
                        end
                    end

                    h(h>1)=1;
                    h(h<0.02)=0.02;

                    %Mean of approximated solutions
                    h = (h_lin+h)/2;

                    %U0=U0+dU/2;dgU0=dgU2;gV0=gV2;dgV0=dgV;
                    %U0=U;dgU0=dgU;gV0=gV;dgV0=dgV;
                    U=U0+dU*h;
                    [gV,dgU] = feval(G_fun,U,varargin{:});
                    if echoReliabPlot
                        plot3(U(1),U(2), gV,'.m')
                    end

                    V = transNf(U',Juz,-ft,Nv,parmeters)';
                    Jzv = inv_normeq(Mu,Su,U,ft);
                    dgV = Jzu*Jzv*dgU;

                end
                %V=V+dV*h;
                %U=V'*Jvu + Mu;

                %% Convergence Test 
                normg=max(norm(gV),abs((dgV'*V)/norm(dgV)/norm(V)-1));
                %verf_grad=norm(dgV/norm(dgV) - V/norm(V));
                if echoReliab 
                    fprintf('Iter %d, Norm G constraint: %g, beta = %g \n', i, max(gV),-norm(V)*sgnG)
                end
                if i>100
                    % STOP
                    normg=0;
                    verf_grad=0;
                end
                i=i+1;
                Us(i+1,:)=U;
                Vs(i+1,:)=V';
                %dUs(i,:)=dgU/norm(dgU);
                %dVs(i,:)=dgV'/norm(dgV);
                gs(i+1)=gV;
            end % end MPP search
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ndgV = norm(dgV);
            if ndgV==0
                beta = -sgnG*(norm(V) - gV);
            else
                beta = -sgnG*(norm(V) - gV/ndgV);
            end

            PF=acu_aprx(-beta);
            if echoReliab
                fprintf('Finish: Iter %d, Norm G constraint: %g, beta = %g \n', i, max(gV),beta)
                figure,plot([Us,gs'],'.-')
                for itex=1:length(U)
                    legtxt{itex}=['U' int2str(itex)];
                end
                legtxt{itex+1}='G';
                legend(legtxt{:});
            end

            %Structural Reliability Methods
            % O. Ditlevsen and H.O. Madsen

            grad_calculation = 1;
            %[~,dgU,varargout] = feval(G_fun,U,varargin{:});
            %dsigdu = varargout{14};

            var_type=1;
            linStressComputed=0;
            [Gs,dgx,varargout] = feval(G_fun,U,varargin{:});
            var_type=0;
            %dsigdx = varargout{14};

            grad_calculation = 0;
            dbeta = -dgx/norm(dgV);
            %P = -cdf(b), dP = -pdf(b)*dbeta;
            dPi = -normpdf(-beta)*dbeta;
        end
        
            
        if i_contin, i_contin=2;end
        
        if ~DinamcVar
            V_pre_optimum = V;
        end
        if echoReliabPlot
            plot3(Us(:,1),Us(:,2),gs','.-r')
            disp([beta, PF, U'])
        end