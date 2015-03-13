function err=vect_err_calc(vr,va,mode)

% Used in the convergence study RBM
% TrecNLOPT.m file

nap=size(va,2);
% nullvr=(vr==0);
% vr(nullvr)=1;
err=zeros(1,nap);
for i=1:nap
%     va(nullvr,i)=va(nullvr,i)+1;
    % Relative error vector
%     ervapi=1-va(:,i)./vr;
    switch mode
    case 1
        err(i)=norm(va(:,i)-vr(:))/norm(vr);
    case 2
        [eri,k]=max(abs(va(:,i)-vr(:)));
        if vr(k)==0; vr(k)=1; end
        err(i)=eri/abs(vr(k));
    end
end