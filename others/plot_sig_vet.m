function plot_sig_vet(area,F,u,sig,tp,bars_out,~)

global X Y Z conect Truss3D ngln
conect0=conect;
if nargin>5
    sig(bars_out)=[];
    conect0(bars_out,:)=[];
end

x=area;
%   Plotagem
%figure
hold on
%fun='tb64';x=[1 1 2];
nb=size(conect0,1);
%axes(findobj(gcf,'Tag','axes4'));

if  Truss3D
    plotfunc = 'plot3';
    xlabel('x');ylabel('y');
    grid on
else
    plotfunc = 'plot';
end    

if tp<0
    % PLOT UNDERFORMADE
    for i=1:nb
        if  Truss3D
            bari = {X(conect0(i,:)),Y(conect0(i,:)),...
                Z(conect0(i,:))};
        else
            bari = {X(conect0(i,:)),Y(conect0(i,:))};
        end
        feval(plotfunc,bari{:},'-.k');
    %plot(X(conect0(i,:)),Y(conect0(i,:)),'-.k');
    end
end
n=length(X);
mx=sum(x)/nb;
esc=4/mx;

if  Truss3D
    glpn = 3;
    escu=max(abs([X(:); Y(:); Z(:)]))/100;
else
    glpn = 2;
    escu=max(abs([X(:); Y(:)]))/100;
end
    escu = 1;
dd=zeros(1,glpn*n);
dd(ngln)=escu*u;

Xd=X(:)+dd(1:glpn:end)';
Yd=Y(:)+dd(2:glpn:end)';
if  Truss3D
    Zd=Z(:)+dd(3:glpn:end)';
end
    
%e=min([sum(abs(diff(X))) sum(abs(diff(Y)))])/2;
% ff=zeros(2*n,1);
% ff(ngln)=F;
% px=ff(1:2:end);
% py=ff(2:2:end);
% ipx=find(px)';
% ipy=find(py)';
%     ex=e/max(abs(px));
%     ey=e/max(abs(py));
% 
% hold on
% for i=ipx
%     quiver (Xd(i),Yd(i),px(i),0,ex/50,'b');
% end
% for i=ipy
%     quiver (Xd(i),Yd(i),0,py(i),ey/50,'b');
% end

medbx = (-X(conect0(:,1))+X(conect0(:,2)))/1.7;
medby = (-Y(conect0(:,1))+Y(conect0(:,2)))/1.7;
if abs(tp)==1
    mat=max(sig);
    mit=min(sig);
    f=sig;
    title('Stress')
    lw=2/mean(x).*x;
elseif abs(tp)==2
    mat=max(sig(:).*x(:));
    mit=min(sig(:).*x(:));
    f=sig(:).*x(:);
    title('Esforços')
    lw=2*x./x;
elseif tp==0
    f=sig;
    title('Stress')
    lw=2/mean(x).*x;
end
if max(lw)>30
    lw=30/max(x).*x;
end
if mat==0, mat=1; end
i=0:30;
ss(i+1,1)=1-i/30;
ss(i+1,2)=i/30;
ss(i+31,2)=1-i/30;
ss(i+31,3)=i/30;
%if tp~=0
%end
for i=1:nb
    %if f(i)>0
    %    c='b';
    %else
    %    c='r';
    %end
    
    if  Truss3D
        bari = {Xd(conect0(i,:)),Yd(conect0(i,:)),...
            Zd(conect0(i,:))};
        % text(medbx(i)+X(conect0(i,1)),Y(conect0(i,1))+medby(i)...
        %,mean(Zd(conect0(i,:))), num2str(sig(i),2));
    else
        bari = {Xd(conect0(i,:)),Yd(conect0(i,:))};
    end
    if f(i)>=0
        %c='b';
        c=[0 1-f(i)/mat f(i)/mat];

    elseif f(i)<0
        %c='r';
        c=[f(i)/mit 1-f(i)/mit 0];

    end
    li=feval(plotfunc,bari{:},...
        'LineWidth',lw(i),'Color',c);
    %text(medbx(i)+X(conect0(i,1)),Y(conect0(i,1))+medby(i)...
    %    ,num2str(sig(i),2));
end

lx=get(gca,'XLim');
mnx=lx(1); mmx=lx(2);
ly=get(gca,'YLim');
mny=ly(1);
mmy=ly(2);
if mmx-mnx>mmy-mny
    mmy=(mmy+mny)/2+(mmx-mnx)/2;
    mny=(mmy+mny)/2-(mmx-mnx)/2;
else
    mmx=(mmx+mnx)/2+(mmy-mny)/2;
    mnx=(mmx+mnx)/2-(mmy-mny)/2;
end
set(gca,'XLim',[mnx mmx]);
set(gca,'YLim',[mny mmy]);
legend
hold off


colorbar('Clim',[mit mat],'YLim',[1 61],'YTick',[1 10 20 30 40 50 61],...
    'YTickLabel',{num2str(mit,4);num2str(2*mit/3,4);...
        num2str(mit/3,4);num2str(0,4);num2str(mat/3,4);...
        num2str(2*mat/3,4);num2str(mat,4)});
%caxis('auto');
colormap(ss);