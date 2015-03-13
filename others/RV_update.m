function updated_varg = RV_update(x,RVkind,RVwich,components,vararg0)
% Update the random variabel vararg0 of the type RVkind 
% with 'components' given as RVwichand using x as new values (or factor).

global link 

    updated_varg = vararg0;
    if RVwich==1||RVwich==51
        % component determinated by design variable
        comps=[];
        for ic=1:length(components)
            comps=[comps; find(link(:,1)==components(ic))];
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Random Variables Type (U)
    switch RVkind
    case 1
        % Design Random Variables
        if RVwich==5
            % Design Variable Factor!!!
            % component determinated by design variable
            comps=[];
            for ic=1:length(components)
                comps=[comps find(link(:,1)==components(ic))];
            end
            if x<1e-9,x=1e-9;end
            updated_varg(comps)=vararg0(comps)*x;
        else
            updated_varg(comps)=x;
        end

    case 2
        % Material Random Variables
        E=vararg0{3};
        if RVwich==1||RVwich==51
            components = comps;
        end
        if RVwich==5||RVwich==51
            E(components)=E(components)*x;
        else
            E(components)=x;
        end
        updated_varg{3}=E;

    case 3
        % Load Random Variables
        fext=vararg0;
        if RVwich==1
            disp('caso n esperado')
        elseif RVwich==2
            disp('caso n esperado')
        elseif (RVwich==3)
            % component determinated by node (using DoF )
            sgnF=sign(components);
            components=abs(components);
            fext(components)=sgnF*x;
        elseif RVwich==4
            % component determinated by DoF 
            sgnF=sign(components);
            components=abs(components);
            fext(components)=fext(components)*0+sgnF(:)*x;
        elseif RVwich==5
            % Factor for all glb, component not needed
            fext=fext*x;
        end
        updated_varg=fext;
    otherwise
            disp('Not recognized RandVarprop type!! (RV_updated)')
        %...others Random Variables
    end
end