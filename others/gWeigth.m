
function b=gWeigth(Xs,xi)
% Comput the weigths of an 
% "gravitational aproximation"
% based on the distance of points precomputed
% Variables:
% Xs - precomputed points
% xi - Desired new point
% b  - weigth coefficients

    % num. of precomputed points
    [np,~]=size(Xs);
    
    % Comput point normalized dist vector 
    dx=norm_diffP(Xs,xi,np);
    
    % Dist norm
    ndx=sum(dx.^2,2).^2;
    
    nz=find(ndx<1e-12);
    if isempty(nz)
        g = 1./ndx;
        Tnd=sum(g);
        b = g/Tnd;
    else
        b = zeros(np,1);
        b(nz(1))=1;
    end
end

function [dx]=norm_diffP(Xs,xi,np)

    % Dist vector Calculation
    dx=Xs-ones(np,1)*xi;
    
    % NORMALIZATION
    % min - max
    xU=max(Xs);
    xL=min(Xs);
    
    % max delta
    dL=xU-xL;
    dL(dL==0)=1;
    
    % Matrization
    dL=ones(np,1)*dL;
    x0=ones(np,1)*xL;
    
    % Dist normal vector
    dx=(dx-x0)./dL;
end