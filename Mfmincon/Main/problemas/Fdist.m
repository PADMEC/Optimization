function F=Fdist(X)

global funinput
Xf=funinput{1};
nf=funinput{2};
for i=1:nf
    F(i)=norm(X-Xf(i,:));
end

