


clear variable global

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CVT Tol
global toll papern
toll=1e-5;
rand('seed',654321)
randn('seed',246810)
%path(pathdef);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set if the gradients are given
global GradOutput_on
GradOutput_on=0;

path('Mfmincon',path) 


path('problemas',path) 

%fun='fun3mymo';%N
%con='con3mymo';
%x0=[1 1 1]*pi/2;

% X
% fun='fun33mo';
% con='con33mo';
% x0=[-1 1 1]*1;

% fun='funRaizmo';
% con='conRaizmo';

% fun='fun3tapmo';
% con='con3tapmo';

% fun='fun30mo';
% con='con30mo';
% x0=[1 1];

% fun='fun3KoskTruss';
% con='con3KoskTruss';
% x0=[2 2 2];

%EXEMPLE 1
% fun='fun_nxmo';
% con='con3xmo';
% x0=[1 1 1]*2;

%EXEMPLE 2
fun='fun33mo';
con='con33mo';
x0=[1 1 1]*1;


%EXEMPLE 3
% fun='Fdist';
% con='ConstrDist';
% x0=[0 0];
% MO_Geom_ProbDef

%EXEMPLE 4
% fun='fun_nxmo';
% con='con4xmo';
% x0=[2 2 2 2]/2;


ptsp=10;
%opts=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Optimizations Types: 'ws','mmx','nbi','nnc','nbim','nncm'
MOptTypes={'ws','mmx','nbi','nnc','nbim','nncm'};

imethods=1:6;
for i=imethods
    MOptype=MOptTypes{i}

    %Run MO method
    [x, f, Cnvrg, fcount,Tot_Time]=Mfmincon(x0, x0*(-inf), x0*inf,fun, con, ptsp, MOptype);

    disp(MOptype)

    %Save method results
    Resum{i}={f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp};
    %[TTs(i),TFCs(i),eveness(i),Par_area1(i),Par_area(i),nnonP(i)]=post_proc(Resum(i){:});
end

% Save Resum = {f,x,Cnvrg,fcount,Tot_Time,MOptype,ptsp}
save ResumOUT Resum

for i=imethods
    [TTs(i),TFCs(i),eveness(i),Par_area1(i),Par_area(i),nnonP(i),ndpP(i)]=post_proc(Resum{i}{:});
    title(MOptTypes{i})
end

%PNG_n_CLOSE
Resume_Tablet=[
TTs;TFCs;nnonP;eveness;Par_area];
%plot Non Domain Graphic
%G_nonD_GRAPHCS

fprintf('Results: method,    Time, F Count, Num. NonPareto, Evness, Pareto Area \n')
for il=imethods
    fprintf('             %d,  %g,     %d,    %d,       %f, %f \n',il,TTs(il),TFCs(il),nnonP(il),eveness(il),Par_area(il))
end

disp(fun)
fprintf('Results: method,    Time, F Count, Num. NonPareto, Evness, Pareto Area \n')
for il=imethods
    fprintf('%d & %.4g & %d & %d & %d & %.4g & %.5g \n',il,TTs(il),TFCs(il),nnonP(il),ndpP(il),eveness(il),Par_area(il))
end

disp([imethods' Resume_Tablet'])
% superf(f);
% figure
% scatter3(f(:,1),f(:,2),f(:,3),30,f(:,3),'filled')
% bf=b*ft;
% hold on,plot3(bf(:,1),bf(:,2),bf(:,3),'x')
% hold off
% figure
% na0=find(Cnvrg>0)+3;
% scatter3(f(na0,1),f(na0,2),f(na0,3),40,f(na0,3),'filled')
% na0=find(Cnvrg==0)+3
% hold on,plot3(f(na0,1),f(na0,2),f(na0,3),'xr')
% scatter3(f(1:3,1),f(1:3,2),f(1:3,3),40,f(1:3,3),'filled')
