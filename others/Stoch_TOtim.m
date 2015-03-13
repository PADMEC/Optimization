function [fre,gre,fob,gob]=Stoch_TOtim(tres,tobj, v_out0, M_out, S_out, dv_out0, dM_out, dS_out)
% [fre,gre,fob,gob]=Stoch_TOtim(tres,tobj, v_out0, M_out, S_out, dv_out0,
% dM_out, dS_out)
%
% Comput Constraint and Objective Values and Gradient
%u,Mu,Su,tmpr,Mtmpr,Stmpr, sn,Msn,Ssn, sig,Msig,Ssig
% u=Mu;
% tmpr=Mtmpr;
% sn=Msn;
% sig=Msig;
% vol=Mvol;

global  mo itp obj_type i_obj

    
    global ResType n_sigma sign_const sign_obj specfc_FF

    fre=[];
    gre=[];
    nval=size(dM_out{1},1);
%ndvab=length(dsn);
%itp=tres;
% ResType = ['probability', 'deterministic']
% sign_const = [1,-1, factor]
% n_sigma(iconst) - R = m+'n'*sig
if isempty(specfc_FF)
    specfc_FF(itp)=0;
end
    for i=itp
        if strcmp(ResType(i),'probability')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if specfc_FF(i)>0
                k=specfc_FF(i);
                % Failure Function specify type
                F = sign_const(i)*(abs(M_out{i}(k)) + ...
                    n_sigma(i)*S_out{i}(k));
                
            % Constraint Gradient
            const_sgn = sign_const(i)*sign(ones(nval,1)*M_out{i}(k)');
            grei = dM_out{i}(:,k).*const_sgn +...
                    n_sigma(i)*dS_out{i}(:,k);
             
            else
                %Trying paper RO transmition truss results
%                 F = sign_const(i)*((M_out{i}) + n_sigma(i)*S_out{i});
%                 const_sgn = sign_const(i)*ones(size(dM_out{i}));

                % Constraint Value
                F = sign_const(i)*(abs(M_out{i}) + n_sigma(i)*S_out{i});

                % Constraint Gradient
                const_sgn = sign_const(i)*sign(ones(nval,1)*M_out{i}');
                grei = dM_out{i}.*const_sgn +...
                        n_sigma(i)*dS_out{i};
            end
            %[frei,cont_max] = max(F);
            frei = F;
               
                
                
            %ibar_test=1:25;
            %ibar_test(13)=[];
            %frei(13) = 0;
            %grei(:,13)=grei(:,13)*1e-6;
            
        elseif strcmp(ResType(i),'deterministic')
            % Constraint Value
            [frei,cont_max] = max(abs(v_out0{i}));
            
            % Constraint Gradient
            const_sgn = sign(v_out0{i}(cont_max));
            grei = dv_out0{i}(:,cont_max)*const_sgn;
            
        end
        fre=[fre;frei];
        gre=[gre grei];
    end
    
            
    if tobj == 0%MULTI-OBJETIVO
        %%%%%%%%    MULTI-OBJETIVO  %%%%%%%%
        tobjs=find(mo(1,:));
    else
         % SCALAR OBJECTIVE
        tobjs=tobj;
    end
    
    nobj=length(tobjs);
    fob=zeros(1,nobj);
    gob=zeros(nval,nobj);
    nM = length(M_out);
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    % Tirar!!!!!!!!!!!!!!!
    % Robust Optimization settings
    %   select same component (max(mean)) 
    %   to optmize the mean and std
    RO_item=zeros(nM,1);
    RO_type=zeros(nM*2,1);
    %for iRO=1:nM
    %    [obsel]=intersect(tobjs,[iRO,iRO+nM]);
    %    if length(obsel)==2
    %        if strcmp(obj_type{iRO},obj_type{iRO+nM})
    %            RO_type(iRO)=1;
    %            RO_type(iRO+nM)=-1;
    %        end
    %    end
    %end
    
    % Begin Objectives Selections
    for i=1:nobj
        % Objective Value and Gradient
        obj = tobjs(i);
        if obj<=nM
            % Mean
            v_obj = M_out{obj};
            dv_obj = dM_out{obj};
        elseif obj>nM
            % STDeviation
            v_obj = S_out{obj-nM};
            dv_obj = dS_out{obj-nM};
        end
        
        if strcmp(obj_type{obj},'specific')
            v_obj = v_obj(i_obj(i));
            dv_obj = dv_obj(:,i_obj(i));
%         else
%             if strcmp(obj_type{obj},'max')
%             %disp('choose an obj_type');
        end
        
        % Minimizatio/Maximization selection "sign_obj"
        % standard option: Minimization
        if isempty(sign_obj)
            signo=1;
        else
            signo=sign_obj(i);
        end
        
        % Robust optimization STD pre-check
        if RO_type(obj)==-1
            cont_max=RO_item(obj-nM);
            v_obj = v_obj(cont_max);
            dv_obj = dv_obj(:,cont_max);
        end

        [fobi,cont_max] = max(signo*abs(v_obj));
        
        % Robust optimization MEAN post-check
        if RO_type(obj)==1
            RO_item(obj)=cont_max;
        end
        gobi = signo*dv_obj(:,cont_max)*sign(v_obj(cont_max));
        
        fob(i)=fobi;
        gob(:,i)=gobi;

        [nx,nv]=size(dv_obj);
        %nv=1;
        if nv>1
        % Constraint others values 
        % G = ((Fj-Fi)^2<0, i -> max value, j -> others)
            
            v_obj(cont_max) = [];
            dv_obj(:,cont_max) = [];

            fcdiff=(signo*abs(v_obj) - (fobi));
            dfcdiff=(signo*ones(nx,1)*sign(v_obj)').*dv_obj - gobi*ones(1,nv-1);
            frei = -fcdiff.^2;
            grei = -(ones(nx,1)*(2.*fcdiff')).*dfcdiff;
            fre=[fre;frei];
            gre=[gre grei];
        end

    end
    
end
        

    