function plotes(arq,tp)

if nargin==1
    tp=1;
end

%   Plotagem
%fun='tb64';x=[1 1 2];
nb=size(arq.barra,1);
%axes(findobj(gcf,'Tag','axes4'));
for i=1:nb
    plot(arq.nox(arq.barra(i,:)),arq.noy(arq.barra(i,:)),'-.k');
    hold on
end
n=size(arq.nox,2);

dd=ones(1,2*n);
dd(arq.ay*2)=0;
dd(arq.ax*2-1)=0;
if ~isfield(arq,'arex')
    arq.arex=ones(1,nb);
end
if ~isfield(arq,'x')
    arq.x=arq.x0;
    %arq.fdis=arq.idis;
    arq.x=arq.x0((arq.arex));
end
if ~isfield(arq,'esc')
    arq.esc=1;
end
dd(find(dd==1))=arq.esc*arq.fdis;
for i=1:n
    arq.noy(i)=arq.noy(i)+dd(i*2);
    arq.nox(i)=arq.nox(i)+dd(i*2-1);
end
%if size(arq.x,2)==nb&max(arq.arex)<nb
%    arq.arex=1:nb;
%end
mx=mean(arq.x);
esc=4/mx;
%medbx = (-arq.nox(arq.barra(:,1))+arq.nox(arq.barra(:,2)))/1.7;
%medby = (-arq.noy(arq.barra(:,1))+arq.noy(arq.barra(:,2)))/1.7;
x=arq.x;
if tp==1
    mat=max(arq.ftes);
    mit=min(arq.ftes);
    f=arq.ftes;
    title('Tensão')
    lw=3/mean(x).*x;
elseif tp==2
    mat=max(arq.ftes(:).*x(:));
    mit=min(arq.ftes(:).*x(:));
    f=arq.ftes(:).*x(:);
    title('Esforços')
    lw=4*x./x;
elseif tp==0
    f=arq.ftes;
    title('Tensão')
end
if max(lw)>30
    lw=30/max(x).*x;
end

i=0:30;
ss(i+1,1)=1-i/30;
ss(i+1,2)=i/30;
ss(i+31,2)=1-i/30;
ss(i+31,3)=i/30;
%if tp~=0
colormap(ss);
colorbar('Clim',[1 61],'YLim',[1 61],'YTick',[1 10 20 30 40 50 61],...
    'YTickLabel',{num2str(mit,4);num2str(2*mit/3,4);...
        num2str(mit/3,4);num2str(0,4);num2str(mat/3,4);...
        num2str(2*mat/3,4);num2str(mat,4)});
%caxis('auto');
%end
for i=1:nb
    %if arq.ites(i)>0
    %    c='b';
    %else
    %    c='r';
    %end
    if f(i)>=0
        %c='b';
        c=[0 1-f(i)/mat f(i)/mat];
            li=plot(arq.nox(arq.barra(i,:)),arq.noy(arq.barra(i,:))...
                ,'LineWidth',lw(i),'Color',c);
            %set(li,);
    elseif f(i)<0
        %c='r';
        c=[f(i)/mit 1-f(i)/mit 0];
            li=plot(arq.nox(arq.barra(i,:)),arq.noy(arq.barra(i,:))...
                ,'LineWidth',lw(i),'Color',c);
            %set(li
    end
    hold on
    %text(medbx(i)+arq.nox(arq.barra(i,1)),arq.noy(arq.barra(i,1))+medby(i)...
    %    ,num2str(arq.ftes(i),2));
end
e=min(sum(abs([diff(arq.nox);diff(arq.noy)])))/2;
ipx=find(arq.px);
ipy=find(arq.py);

for i=ipx
    ex=e/max(abs(arq.px));
    quiver (arq.nox(i),arq.noy(i),arq.px(i),0,ex,'b');
end
for i=ipy
    ey=e/max(abs(arq.py));
    quiver (arq.nox(i),arq.noy(i),0,arq.py(i),ey,'b');
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