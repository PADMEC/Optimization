

L=3;
Xs=[-1,-1;-1,1;1,-1;1,1; lhsdesign(12,2)*2-1]*L;
[np,ndv]=size(Xs);
Ys=zeros(np,1);
for i=1:np
    Ys(i,:)=peaks(Xs(i,1),Xs(i,2));
end

xi=[.16,.32]*L;

[yi1]=aprox1(Xs,Ys,xi);
[yi2]=aproxg(Xs,Ys,xi);
yi3 = griddatan(Xs,Ys,xi);

Yv=peaks(xi(:,1),xi(:,2));

disp([Yv, yi1, yi2,yi3])

xi=(-1:0.1:1)*L;
[X1,X2]=meshgrid(xi,xi);
X=[X1(:),X2(:)];
%b=Xs'\X';
%Ya=Ys'*b;
Yv=peaks(X(:,1),X(:,2));

nx=length(X);
Ya1=zeros(nx,1);Ya2=Ya1;
for i=1:nx
    xi=X(i,:);
    %[yi]=aprox1(Xs,Ys,xi);
    %Ya1(i,:)=yi;
    
    [yi]=aproxg(Xs,Ys,xi);
    Ya2(i,:)=yi;
    
    %yi = interpn( Xs(:,1), Xs(:,2),Ys, xi(1),xi(2));
    yi = griddatan(Xs,Ys,xi);
    %interp1(x,Y,xi)
    Ya3(i,:)=yi;
end
%[norm(Yv-Ya1),norm(Yv-Ya2),norm(Yv-Ya3)]
%[max(abs(Yv-Ya1)),max(abs(Yv-Ya2)),max(abs(Yv-Ya3))]
[norm(Yv-Ya2),norm(Yv-Ya3)]
[max(abs(Yv-Ya2)),max(abs(Yv-Ya3))]


figure
grid on
Yvm=peaks(X1,X2);
mesh(X1, X2, Yvm)
hold on, hidden off
plot3(Xs(:,1), Xs(:,2), Ys,'o')
plot3(X(:,1), X(:,2), Ya2','.m')

figure
Yvm=peaks(X1,X2);
mesh(X1, X2, Yvm)
hold on, hidden off
plot3(Xs(:,1), Xs(:,2), Ys,'o')
plot3(X(:,1), X(:,2), Ya3','.k')


Yvm=peaks(X(:,1), X(:,2));
err2=abs(Ya2'-Yvm');
err3=abs(Ya3'-Yvm');

figure, hold on,
plot(err2-err3,'.-r')

figure, hold on,
plot(abs(Ya2'-Yvm'),'.-r')
plot(abs(Ya3'-Yvm'),'.-k')

figure, hold on,
hold on, hidden off
plot3(X(:,1), X(:,2), err2,'.r')
plot3(X(:,1), X(:,2), err3,'.k')
