function  [fn,kni,Z,kvi,Z0,kdi,kmni,Zd]=feoffline2(SN,tipo)
%
%
%  offline stage for FE analysis:
%
%
%   calculates the submatrices and load vector 
%
%
%  solves for each dv in turn

load trussdata.mat
N=size(SN,1);

nmf=tipo(2);
nmv=tipo(3);
kvi=[];Z0=zeros(0,0,0);
kdi=[];kmni=[];Zd=zeros(0,0,0);
for in = 1:N
%LargTruss
        %mu = [SN(in,:) 1]
        mu = SN(in,:);
        K=K1{1}*0;
        for i=1:size(K1,2)
            K = K+mu(i)*K1{i};
        end

    if tipo(1)>0 |tipo(2)>0

    %
    %   get solutions for each dv pair 
    %
        %x  = Ksolve2(mu);
        x = full(K\F);

    %   build the vector of all solutions
    %
        Z(:,in)= x;
    end
    
    if tipo(2)>0
        
        %%Linear Buckling
        
        area=mu(link(:,1));
        esf = tresf(area,comp,els,ang,glb,x);
        k0 = kgeom(esf,comp,ang,glb);
        
        [V,D]=eig(full(K),full(k0));
        
        D=diag(D);
        [Lambd,mdf]=sort(abs(real(D)));
        du=full(V(:,1:mdf));
        Z0(:,in,1:nmf)=du(:,1:nmf);
    end
    
    if tipo(3)>0
        
        %Vibration Analysis
        
        area=mu(link(:,1));
        Km=Km1{1}*0;
        for i=1:size(Km1,2)
            Km = Km+mu(i)*Km1{i};
        end
        
        [D,W2]=eig(full(K),full(Km));
        
        [Lambd,mdv]=sort(abs(real(W2)));
        
        du=full((D(:,1:mdv)));
        Zd(:,in,1:nmv)=du(:,1:nmv);

    end
 %
 %  K*Z matrices
 %
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPACE CHECK
%lld=lindep(Z');
lld=[];
if ~isempty(lld)
    disp('Z linearmente dependente, colunas:')
    disp(num2str(lld))
    Z(:,lld(:,2))=[];
    
    disp('SN-considerado')
    SN(lld(:,2),:)=[]
    %warning(num2str(lld))
end

lld=lindep([]);%Z0');
if ~isempty(lld)
    disp('Z0 linearmente dependente, colunas:')
    disp(num2str(lld))
    %warning(num2str(lld))
end

lld=lindep([]);%Zd');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(lld)
    disp('Zd linearmente dependente, colunas:')
    disp(num2str(lld))
    %warning(num2str(lld))
end


 %
 %  Kn matrices
 %
for i=1:size(mu,2)
   kni(:,:,i) = full(Z'*K1{i}*Z);
end

 %
 % the vector fn
 %
    fn  = Z'*F;

for imdf=1:size(Z0,3);
    z0i=Z0(:,:,imdf);
    %k0n(:,:,imdf)=z0i'*k0*z0i;
    for idv=1:size(mu,2)
       kvi{idv}(:,:,imdf) = full(z0i'*K1{idv}*z0i);
    end

end

for imdv=1:size(Zd,3);
    zdi=Zd(:,:,imdv);
    for idv=1:size(mu,2)
       kdi{idv}(:,:,imdv) = full(zdi'*K1{idv}*zdi);
    end
    for idv=1:size(mu,2)
       kmni{idv}(:,:,imdv) = full(zdi'*Km1{idv}*zdi);
    end
end

        %kvi=>{n da variavel}(gral de liberdade,n da amostra, n do modo d falmbagem)
