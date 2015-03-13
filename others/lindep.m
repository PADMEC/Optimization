function [ldep,jcontd]=lindep(A)

[nl,nc]=size(A);
%if nc>nl
%    A=A';
%    [nl,nc]=size(A);
%end

nA=null(A');
[nl,nnul]=size(nA);
id_out=[];
for i=1:nnul
    [sA,indx]=sort(abs(nA(:,i)),'descend');
    [intrs, i_id]=intersect(indx,id_out);
    ids_ok=1:nl;
    ids_ok(i_id)=[];
    id_out=[id_out indx(ids_ok(1))];
end

ldep=[id_out' id_out'];
jcontd=[id_out];
A=full(A);
for i=1:nl
    if sum(find(jcontd==i))==0
        for j=i+1:nl
            if sum(find(jcontd==j))==0
                %tet=subspace(A(i,:)',A(j,:)')/norm((A([i,j],:)));
                kmed=mean(A(i,:))/mean(A(j,:));
                tet=norm( A(j,:)*kmed-A(i,:) )/norm(A(i,:));
                if tet<1e-10
                    ldep=[ldep;i j];
                    if sum(find(jcontd==i))==0
                        jcontd=[jcontd i j];
                    else
                        jcontd=[jcontd j];
                    end
                end
            end
        end
    end
end