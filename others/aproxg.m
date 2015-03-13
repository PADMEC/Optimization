
function [yi]=aproxg(Xs,Ys,xi)

    b=gWeigth(Xs,xi);
    yi=Ys'*b;
end
