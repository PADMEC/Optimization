
global G_const G_limits betatarg RBtype linStressComputed uniq_Failure
global Stochast_on RandVarProp rv_data_comp RandVarDist Correl DinamcVar Method
global serieLimitState Stochast_moment var_type GradOutput_on
linStressComputed=0;
GradOutput_on=1;

 load  NewPortcLarg2F
 %arq.x0=[0.1841,1.0218,1.3823];%Lin FEM
 %arq.x0=[ 0.1820,1.0136,1.3914];%Lin POD
 %arq.x0=[ 0.2104,1.0396,1.3516];%Lin POD
 
 arq.x0=[1 1 1];
 %arq.x0=[0.2871    0.9702    1.3948];
 %arq.x0=[0.1990    0.9764    2.3947];
 %0.1769    0.9814    1.3134];
 %        0.1769    0.9814    1.3140
 
 % RBDO optimum
%arq.x0 = [0.2188    0.9725    2.3901];
%arq.x0 = [0.2192    0.9735    2.3888];

%FEM RBDO
% NEW
%arq.x0 = [0.2987, 1.0128, 2.4630];
% ffobj = 1.4541;
% iterations: 11
% funcCount: 26
% algorithm: 'sequential quadratic programming'
% message: [1x782 char]
% constrviolation: 1.2185e-010
% stepsize: 1
% firstorderopt: 1.3867e-007
% tempo de otimizacao:  266.9341

%POD RBDO
% NEW
% arq.x0 = [0.3007, 1.0103, 2.4625];
% ffobj = 1.4534;
% iterations: 11
% funcCount: 26
% algorithm: 'sequential quadratic programming'
% message: [1x782 char]
% constrviolation: 0
% stepsize: 1
% firstorderopt: 1.4414e-007
% tempo de otimizacao:   78.2626
%POD Dual procedure
%arq.x0 = [0.2066    0.9944    2.4209];

%POD R2DO
%  arq.x0 =[ 0.7326    1.4872    3.1350]/2;
%  arq.x0 =[ 0.4096    1.3092    3.4711]/2;
      arq.x0 = [0.4595    0.4029    1.8587];%Mean
%      arq.x0 = [0.3833    0.7498    1.5529];%STD
%     arq.x0 = [0.4189    0.5845    1.6990];%average

   X_ops = arq.x0;

%MPa = 1MN/(1e2cm)^2 = 1MN/1e4cm^2 = 1e2N/cm^2
%2e5MPa = 2e7 N/cm^2
arq.xmin=0.1;
arq.xmax=10;
arq.E=2.07e7;
E=arq.E;dE=E/5;% 5,0%
%k_p =6000;
k_p = 2000;%500
RBtype=1;N=60;
% RBtype=0;

Stochast_on=1;
%DinamcVar = 1e5;Method='MC';
%DinamcVar = 0;Method='FORM';
% DinamcVar = 0;Method='FERUM';

Stochast_moment = 1;
% DinamcVar = 100000;Method='MC';
% DinamcVar = 3;Method='PC';
DinamcVar = 3;Method='R2DO';




%% Objective
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5 - Lambd
%Inicializing...
n_tot_obj = 5;
global mo  obj_type i_obj
mo = zeros(1,n_tot_obj*2);
obj_type=cell(1,n_tot_obj*2);

%arq.tobj=2;%Displacement
%arq.tobj=3;%stress

iobj=3; %Stress
objs=[iobj,n_tot_obj+iobj];
arq.tobj=0;%Multi_Obj
mo(objs)=1;
 obj_type{objs(1)}='max';
 obj_type{objs(2)}='max';
 i_obj = 2020; % i_obj = 1281;% i_obj = 1241;

% iobj=2;% Displacement
% objs=[iobj,n_tot_obj+iobj];
% arq.tobj=0;%Multi_Obj
% mo(objs)=1;
%  obj_type{objs(1)}='specific';
%  obj_type{objs(2)}='specific';
%  i_obj = 661;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Constraint
global set_Failure_FORM
G_fun = 'StructLimitState';
S_fun = 'Stoch_Analysis';

% Reliability Based Constraint
%G_const = v_out = {vol,u,vonmis,sn,Lambd,dvol,du',dsig',den',dLambd};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5 - Lambd

%Failure Function
 G_const = 2;
 G_limits = 3.5;
 set_Failure_FORM = 661; % top left node - x direction
 
betatarg = 3.3;
serieLimitState = [];

uniq_Failure = 1;


% Others Constraints
global ResType n_sigma sign_const
% 
arq.tpres=[0 0 0 0 0];
iconst = 1; %Vol
arq.tpres(iconst)=1;
ResType{iconst}='deterministic';
n_sigma(iconst)=0;
sign_const(iconst)=1;


arq=MatProp(arq);
% arq.MatCurv(:,2)=arq.MatCurv(:,2)/1e4;
arq.E=arq.MatCurv(1,2)/arq.MatCurv(1,1);
SIGu=arq.MatCurv(end)*.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor

% RandVarDist(1,:)=[2 2 2 1]; % distribuction
% RandVarDist(2,:)=[E E E 1]; % Mu par1
% RandVarDist(3,:)=[dE dE dE .25]; % std par2

%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)


RandVarProp(1,:)=[3  3]; % type var
RandVarProp(2,:)=[4 4]; % which type data

ipx=[1,22,43,64,85,106,127,148,169,190,211,232,...
    243,254,265,276,287,298,309,320,331];
ipy=[331,332,333,334,335,336,337,338,339,340,341];
rv_data_comp = {-ipy'*2,ipx'*2-1};

RandVarDist(1,:)=[1  2]; % distribuction
RandVarDist(2,:)=[k_p*2  k_p]; % Mu  par1
RandVarDist(3,:)=[k_p  k_p/2]; % std par2


%Mu = [1,1,1]*.5;
%Cu = diag([1 2 3]/50);
%ft=[2,2,2];


nrv=size(RandVarProp,2);
Correl=eye(nrv);