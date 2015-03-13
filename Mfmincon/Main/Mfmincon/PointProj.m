
    function [Pp,xvr,tp]=PointProj(Vr,P0,ps)
    % Project 'ps' points in the plane form by 'Vr' vectors
    % and 'P0' is the reference point

    % Pp - are the projected points
    % xvr - are the projected points in the plane coordenates (p=xvr*Vr)
    % tp - are the distances from the proj plane to the originals points

    [ndiv,nfun,npar]=size(ps);
    
    tp=zeros(npar,ndiv);
    Pp=zeros(ndiv,nfun,npar);
    xvr=zeros(nfun-1,ndiv,npar);
    
    %NORMAL VECTOR GENERATE
    Pt=null(Vr');
    %%%%%%%%%%%
    %PROJECTION
    %P0+Vr*abp+Pt*t=pf
    %[Vr Pt]*abt=pf-P0
    %A*abt=B
        A=[Vr -Pt];
        Am1=A^-1;
        for i=1:npar
            B=(ps(:,:,i)'-P0*ones(1,ndiv));
            abt=Am1*B;

            tp(i,:)=abt(end,:);
            Pp(:,:,i)=ps(:,:,i)+(Pt*tp(i,:))';
            xvr(:,:,i)=(abt(1:(end-1),:));
            %pfa(i,:)=(abt(1,:));
            %pfb(i,:)=(abt(2,:));

        end
    end
    