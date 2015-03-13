function R=con3KimB(x)

R= - cos(x(:,1)) - exp(-x(:,2)) + x(:,3);
R=[R -x(:,1:2)];
R=[R 1.2-x(:,3)];
R=[R x(:,1)-pi];