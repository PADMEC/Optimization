function [varargout]=FunStatistcalc(S_fun, varargin)

%     S_fun = 'StructLimitState';
%     S_fun = 'Stoch_Analysis';

    global Stochast_on RandVarDist Correl DinamcVar Method betatarg Stochast_moment
    if Stochast_on
        ft = RandVarDist(1,:);
        Mu = RandVarDist(2,:);
        Su = RandVarDist(3,:);
        Cov= Correl.*(Su'*Su);
        %Su=diag(Cov)'.^.5;
        %C=Cov./(Su'*Su);
        
        N1 = DinamcVar;
%         if strcmp(Method,'FORM')
%                 DinamcVar=i_contin;
%         %else
%         end
        
        if Stochast_moment
        
            [ v_out0, M_out, S_out, dv_out0, dM_out, dS_out ] = StochastiCalc( Mu,Cov,ft,S_fun,Method,N1,varargin{:});
            varargout = {v_out0, M_out, S_out, dv_out0, dM_out, dS_out};
            
        else

            [PF,beta, dbeta] = Reliab(Mu,Cov,ft,S_fun,Method,DinamcVar,1,varargin{:});
            varargout = {PF,beta, dbeta};
            %disp([xvalu, dbeta'])
        end
    end