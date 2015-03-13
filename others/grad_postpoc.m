function [dsig,dF,den,desf]=grad_postpoc(du, u, sig, fext, E, area, link, Tu_strain, Mesf_F, Var_kinds)


    
    %if var_type==2
        % d2f/dudx, Random Variables nd desg variables
        %d(Fi=Fe)/dx,  + (dFi/du)*(du/dx) = dFe/dx
        %dd(Fi=Fe)/dxdy, (ddFi/dxdy) + (ddFi/du2)*(du/dx)*(du/dy)
        % + (dFi/du)*(ddu/dxdy) = dFe/dxdy
        % = K*ddu/dxdy = dFe/dxdy - ddFi/dxdy - dux'*K2*duy
        %F/A = E*u/L
        %K=EA/l*[], dK/du = K2 = 0
        
        %case x=area, y=Fei
        %.:  K*ddu/dxdy = 0 - 0 
        
    %end
    
    ndv=size(du,2);
    % Sensibility outputs
    des = Tu_strain*du;
%     [sig,E]=strain2stress(mat_curv,es);
    dsig = E*ones(1,ndv).*des;
    da= zeros(length(area),ndv);
    for ibv=1:ndv
        if Var_kinds(ibv)==1
            da(link(:,1)==ibv,ibv)=1;
        end
    end
    desf = dsig.*(area'*ones(1,ndv)) + (sig*ones(1,ndv)).*da;
    dF=Mesf_F*desf;
    den=fext'*du ;%+ u'*dF;
end