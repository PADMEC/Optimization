function [yi]=aprox1(Xs,Ys,xi)

    [np,ndv]=size(Xs);

    dx=Xs-ones(np,1)*xi;
    ndx=sum(dx.^2,2);
    
    nz=find(ndx<1e-12);
    if isempty(nz)
        [d1,k1]=min(abs(dx));
        [d,k]=sort(ndx);
        p_sel=[k(1:ndv) k1'];
        xsel=Xs(p_sel,:);
        ysel=Ys(p_sel,:);
        b=xsel'\xi';
        yi=ysel'*b;
    else
        yi=Ys(nz(1),:);
    end
end
