function [x, ft, Cnvrg, fcount,Tot_Time]=Mfmincon(x0,vlb,vub, fun, con, ptsp, MOptype,varargin)
%Optimizations Types: 'ws','mmx','nbi','nnc','nbim','nncm'
%MOptype='NBIm'

MOptTypes={'ws','mmx','nbi','nnc','nbim','nncm'};

method=find(strcmpi(MOptype,MOptTypes));
if isempty(method)
    disp('Mult-obj Opt method not found')
end

frestMetodos={'wscon','mmxcon','nbicon','nnccon','nbicon','nnccon'};
fobjMetodos={'wsfun','mmxfun','nbifun','nncfun','nbifun','nncfun'};
funmtd=fobjMetodos{method};
conmtd=frestMetodos{method};

fcount=[];
Cnvrg=[];

global fs GradOutput_on AVeq BVeq
fs0=fs;
fs=6;
f0=feval(fun,x0,varargin{:});
%fs=fs0;
nfun=length(f0);

    %'TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',10000,...
options=optimset('Display','iter',...
    'LargeScale','off','Algorithm', 'active-set');%'sqp');%'interior-point');%
    %'LargeScale','off','Algorithm', 'interior-point');

  options=optimset(options,'TolX',1e-6,'MaxIter',50,'MaxFunEvals',250,'TolFun',1e-6);    
%     options=optimset(options,'LineSearchType','cubicpoly');
%    options=optimset(options,'LargeScale','on');
%    options=optimset(options,'GradObj','on','GradConstr','on','MaxFunEvals',40000,'Algorithm' ,'active-set');
%    options=optimset(options,'MaxIter',10000,'Display','iter','DerivativeCheck','off', 'LargeScale','on');%,'MaxTime',9);

if GradOutput_on
    options=optimset(options,'GradObj','on','GradConstr','on','DerivativeCheck','off');
end

AUXILIAR_OPTIMIZATION=1;

%%%%%%%%%%%%%%%%%%%%%%%%%
% % HYPLAS TEST
% X_ops = [
%   357.6893  276.1824
%   391.3865  100.0000
%   100.0000  100.0000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0=[1,2,2.5];
 x01=x0;
%  x0 = [0.4695    0.4093    1.8473]+.2;
for i=1:nfun
    % HYPLAS TEST
%      x0=X_ops(i,:);
    [xi,ffobj,Converg,Resumo]=fmincon(funmtd,x0,[],[],AVeq,BVeq,...
        vlb,vub,conmtd,options,fun,con,i,[],[],[],[],[],varargin{:});
    x(i,:)=xi;
    ft(i,:)=feval(fun,xi,varargin{:});
    inv=1:nfun; inv(i)=[];
    %for j=inv
    fprintf('Try scalar opt F%d = %s\n',i,num2str(ft(i,:),4));
    fprintf('Design try: %s\n',num2str(x(i,:)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUXILIAR OPTIMIZATION (Check redundance variables x objective)
    if AUXILIAR_OPTIMIZATION
%         Algopt =optimget(options,'Algorithm');% Get Algorithm
%       options=optimset(options,'Algorithm','sqp');%'interior-point');%'active-set');%
        x0i=xi;
%         options=optimset(options,'TolX',1e-4);
        [xi,ffobj,Convrg,Resm]=fmincon('constfun',x0i,[],[],AVeq,BVeq,...
           vlb,vub,'constcon',options,fun,con,inv,ft(i,i),0,ft(i,:),i,varargin{:});    
%         options=optimset(options,'TolX',[]);

        x(i,:)=xi;
        ft(i,:)=feval(fun,xi,varargin{:});
%         options=optimset(options,'Algorithm', Algopt);%'sqp');%'interior-point');%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       %x0i=xi;
    %end
    %x0i=x0;
%     x0=x0i;
    funcCount=Resumo.funcCount;
    fcount=[fcount funcCount];
    Cnvrg=[Cnvrg Converg];

    fprintf('Final scalar opt F%d =   %s\n',i,num2str(ft(i,:),4));
    fprintf('Convergencia : %d   funcCount: %d \n\n'...
         ,Converg,funcCount);
    fprintf('Final Design : %s\n',num2str(x(i,:)));

    % 2D Large Truss (2 objs)
%         x0=[0.2467    0.7655    1.6050];

end
 x0=x01;
 x0i=x01;
Anchor_points=ft

if 1==0
    %Restriction View
    minx=min(x);
    mxx=max(x);
    lx=mxx-minx+2;
    nd=19;
    xx=minx(1)-1:lx(1)/nd:mxx(1)+1;
    yy=minx(3)-1:lx(3)/nd:mxx(3)+1;
    R=[];
    Ff=[];
    for xxi=xx
        j=0;
        R=[R;zeros(1,nd+1)];
        Ff=[Ff;zeros(1,nd+1,nfun)];
        for xj=yy
            j=j+1;
            Rj=feval(con,[xxi xxi xj]);
            Fj=feval(fun,[xxi xxi xj]);
            R(end,j)=Rj(2);
            Ff(end,j,:)=Fj;
        end
    end
    [X,Y]=meshgrid(xx,yy);
    figure
    contour3(X,Y,R,[-1 0 1])
end

b=betaweigs(nfun,ptsp,1);

par=1:nfun;

[T,mf,l,n,t]=MO_precomputs(ft,method);

%[x,ff] = fgoalattain('nbifun',x0,mf,abs(mf),[],[],[],[],vlb,vub,'nbicon',options,fun,con,0,0,0,0)

f=zeros(size(b));
    
%  options=optimset(options,'Algorithm','active-set');%'sqp');%'interior-point');%

for i=1:size(b,1), bt(i,:)=b(i,:)*ft;end

t0=cputime;
% x0i=[120,120];
if method>4
%     x=x(par,:);
%     f=f(par,:);
    MyMtdfs
else
    diffx=[0,0,0];
    aT=ft;
    %%%%%%%%%%%%
    %xx=0:0.05:2;figure,hold on
    for i=1:size(b,1)
        a=T*b(i,:)';
    %    a=sin(pi/2*a);
    %    ff=a(1)*(2-xx).^.5+a(2)*xx;plot(xx,ff)
        %bt(i,:)=a'.*l+mf;
        if (method==4)||(method==2)
            ax=a.*l'+mf';
        elseif method==3
            ax=a+mf';
        else
            ax=a;
        end
        [mnd,ix]=min(sum((ax*ones(1,size(aT,1))-aT').^2));
        x0=x(end,:);
        % Choose nerest solution
        if norm(ax'-aT(end,:))>1.2*norm(ax'-aT(ix,:))
            x0=x(ix,:);
        end
%         if method==3
%             t=0;
%         end
%         
%           if method==2
%             %Min-Max try
%                x0=x0i;
%           end
        %[AVeq t*0]
        [xt,ffobj,Converg,Resumo]=fmincon(funmtd,[x0 t],[],[],[],BVeq,...
            vlb,vub,conmtd,options,fun,con,-1,mf,l,a,n',par,varargin{:});
        funcCount=Resumo.funcCount;
        fcount=[fcount funcCount];
        Cnvrg=[Cnvrg Converg];
        if method==3||(method==5)
            t=xt(end);
            ts(i)=t;
            xt(end)=[];
        end
        diffx(i)=norm(xt-x0);
        aT(end+1,:)=ax;
        x(end+1,:)=xt;
        f(i,:)=feval(fun,x(end,:),varargin{:});
        
        % Auxiliar optimizations
        if method==1||(method==2)
        if ~isempty(find(b(i,:)==0,1))
            inv=find(b(i,:)==0);
            par=find(b(i,:));
            options=optimset(options,'TolCon',1e-4);
            [xi,ffobj,Convrg,Resm]=fmincon('constfun',xt,[],[],AVeq,BVeq,...
             vlb,vub,'constcon',options,fun,con,inv,f(i,par),0,f(i,:),par,varargin{:});    
            options=optimset(options,'TolCon',1e-6);
            x(end,:)=xi;
            f(i,:)=feval(fun,x(end,:),varargin{:});
        end
        end
        
        ft(end+1,:)=f(i,:);
        fprintf('F_%d =   %s\n',i+nfun,num2str(f(i,:),4));
        fprintf('Convergencia: %d   funcCount: %d \n\n'...
             ,Converg,funcCount);

    end
end
Tot_Time=cputime-t0;
%f=[ft(1:3,:);f];
%figure,plot(diffx,'.-')

%feval(con,x,'p');