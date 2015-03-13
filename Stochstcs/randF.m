function [z,dz,a]=randF(x,a,t_ord)
global echo
if isempty(x)
    z=[];dz=[];
    return
end

if nargin<3
    t_ord = 0;
end
if t_ord
    u = a;
else
    u = x;
    x = a;
end

d = x;E = 64/(pi^2);
    
plot_restrain=-1;
if plot_restrain==1
    xs=.1:.2:8;
    us=.1:.2:8;
    [x,u]=meshgrid(xs,us);

    %1/2^.5-.5
    
    
    %Pc = pi.E.I/L^2, I = pi.(D^4 - d^4)/64;
    D = d+2*u;L = D.*x;
    I = pi/64*(D.^4 - d.^4);
    z = -pi*E*I./(L.^2) + 10;
    %z = (D.^4 - d.^4)./(L.^2);

    figure,hold on,contour(u,x,z,[0 0])
    dz=0;
    return
elseif (echo)&&(plot_restrain==0)
        plot(u(1),u(2),'.')
end
D = d+2*u(:,1);L = D.*u(:,2);
I = pi/64*(D.^4 - d.^4);
z = -pi*E*I./(L.^2) + 10;

%z2 = (D.^4 - d.^4)./(L.^2);

global RandVarProp rv_data_comp var_type
sens = 1;

if sens
%% Sensibility analysis
%   DF/Dx = 0 = dFi0/du * du/dx + dF/dx 

%   syms ds u1 u2 xs
%   Ds = xs+2*u1;
%   Ls=Ds*u2;
%   f=(Ds.^4 - xs.^4)./(Ls.^2);

    if var_type==1
        % Design Variables
        %simplify((diff(f,xs))) = 
        %(8*u1)/u2^2 - (32*u1^4 + 32*xs*u1^3)/(u2^2*Ds^3)
        
        dz = (8*u(1))/u(2)^2 - (32*u(1)^4 + 32*x*u(1)^3)/(u(2)^2*D^3);
    
    else
        % Random Variables nd d2f/dudx
        
%simplify((diff(f,u1))) = (4*Ds)/u2^2 + (4*xs^4)/(u2^2*Ds^3)
        dz(1,1) = -4*D./u(2).^2 - 4*d.^4./(u(2).^2.*D.^3);
        
%simplify((diff(f,u2))) = (2*(ds^4 - (Ds)^4))/(xs^3*(Ds)^2)
        dz(2,1) = -2*(d.^4 - D.^4)./(u(2).^3.*D.^2);
    end
else
    dz=0;
end