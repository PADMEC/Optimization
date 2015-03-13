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


global fs GradOutput_on
fs0=fs;
fs=6;
f0=feval(fun,x0,varargin{:});
%fs=fs0;
nfun=length(f0);

options=optimset('Display','iter','LargeScale','off',...
    'Algorithm','sqp','TolFun',1e-7,'MaxFunEvals',600);

if GradOutput_on
    options=optimset(options,'GradObj','on','GradConstr','on','DerivativeCheck','off');
end

x0i=x0;
for i=1:nfun
    [xi,ffobj,Converg,Resumo]=fmincon(funmtd,x0,[],[],[],[],...
        vlb,vub,conmtd,options,fun,con,i,[],[],[],[],[],varargin{:});
    x(i,:)=xi;
    ft(i,:)=feval(fun,xi,varargin{:});
    inv=1:nfun; inv(i)=[];
    %for j=inv
    
%     options=optimset(options,'TolCon',1e-4);
%        [xi,ffobj,Convrg,Resm]=fmincon('constfun',xi,[],[],[],[]...
%            ,vlb,vub,'constcon',options,fun,con,inv,ft(i,i),0,ft(i,:),i,varargin{:});    
%     options=optimset(options,'TolCon',1e-6);

       disp([x(i,:);xi])
       
%        x(i,:)=xi;
%        ft(i,:)=feval(fun,xi,varargin{:});

       %x0i=xi;
    %end
    %x0i=x0;
    
    funcCount=Resumo.funcCount;
    fcount=[fcount funcCount];
    Cnvrg=[Cnvrg Converg];
end

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
    
%T=eye(nfun);
%n=zeros(1,nfun);
par=1:nfun;
%mf=min(ft);
%l=max(ft)-mf;
%t=[];
%ini=3;

    [T,mf,l,n,t]=MO_precomputs(ft,method);

%options=optimset('Display','iter','LargeScale','off','MaxFunEvals',10000,'MaxIter',100);

%[x,ff] = fgoalattain('nbifun',x0,mf,abs(mf),[],[],[],[],vlb,vub,'nbicon',options,fun,con,0,0,0,0)

%x=[x1;x2;x3];
f=zeros(size(b));

for i=1:size(b,1), bt(i,:)=b(i,:)*ft;end

t0=cputime;

if method>4
    
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
        if norm(ax'-aT(end,:))>2*norm(ax'-aT(ix,:))
            x0=x(ix,:);
        end 
        
        if method==2
            %Min-Max try
            x0=x0i;
        end
        
        [xt,ffobj,Converg,Resumo]=fmincon(funmtd,[x0 t],[],[],[],[],...
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
            [xi,ffobj,Convrg,Resm]=fmincon('constfun',xt,[],[],[],[],...
             vlb,vub,'constcon',options,fun,con,inv,f(i,par),0,f(i,:),par,varargin{:});    
            options=optimset(options,'TolCon',1e-6);
            x(end,:)=xi;
            f(i,:)=feval(fun,x(end,:),varargin{:});            
        end
        end
        
        ft(end+1,:)=f(i,:);
        fprintf('F%d =   %s\n',i,num2str(f(i,:),4));
        fprintf('Convergencia: %d   funcCount: %d \n\n'...
             ,Converg,funcCount);

    end
end
Tot_Time=cputime-t0;
%f=[ft(1:3,:);f];
%figure,plot(diffx,'.-')

%feval(con,x,'p');