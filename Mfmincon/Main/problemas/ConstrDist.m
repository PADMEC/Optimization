function R=ConstrDist(X)

global coninput
Xc=coninput{1};
d0=coninput{2};
nc=coninput{3};

for i=1:nc
    R(i)=d0-norm(X-Xc(i,:));
    
end

