function f=fun33mo(x)

    f(1)=x(1)^3 + x(2) + 2*x(3);
    f(2)=x(2)^3 + x(1) + 2*x(3);
    %f(2)=x(1)^2*x(2)+x(1)^2*x(3);
    f(3)=-prod(x)/125;%x(3)*(x(2)+0)*(x(1)+0);
    
    
%f=f([2 3 1]);