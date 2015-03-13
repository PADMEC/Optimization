%MY CVT
clear all
global toll
toll=5e-4;

n=31;
p=2;
%nt=20;nm=2;ni=0;
ntest=100;

X=amostraRB(p,0,1,n);
%X=rand(n,p);

seed=123456;
sample_function_init=3;
sample_function_cvt=1;
sample_num_cvt=10000;maxit=30;
generator=X';
tic
[ generator_new, seed_new ] = cvt ( p, n, sample_function_init, ...
  sample_function_cvt, sample_num_cvt, maxit, seed, generator );
toc
figure,voronoi(generator_new(1,:),generator_new(2,:));

reset = 1;

%     [ X, seed_new ] = region_sampler ( m, n, n, ...
%       sample_function_init, reset, seed );

X=uniformCVT(X,p,n,100*n);

% figure,hold on
% voronoi(X(:,1),X(:,2));
% [V,C] = voronoin(X);
% for in=1:n
%    ci=C{in};
%    ci(find(ci==1))=[];
%    plot(mean(V(ci,1)),mean(V(ci,2)),'.g')
% end
% plotnum(X(:,1),X(:,2),'.')
% xlim([0 1]),ylim([0 1])
