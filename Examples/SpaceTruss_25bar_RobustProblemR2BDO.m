
global G_const G_limits betatarg RBtype linStressComputed linear_prob
global Stochast_on RandVarProp rv_data_comp RandVarDist Correl DinamcVar Method
global serieLimitState Stochast_moment GradOutput_on
linStressComputed=0;
GradOutput_on=1;

SpaceTruss_25
 
arq.x0=[1 1 1 1 1 1]*2.3;

% %Article Min mean(c)
%  arq.x0=[0.0500	0.0500	5.74	1.718	1.054	5.574];
% %Article Min std(c)
%  arq.x0=[0.147	0.672	3.465	0.566	0.822	8.048];

  %arq.x0=X_ops;
%   arq.x0=X_ops(1,:);

X_ops = arq.x0;


 %arq.x0=[.05 .05 5.74 1.72 1.05 5.57];
 %SD opt [6,9:14][ 0.3788    0.1717    3.8639    0.0526    1.7117    6.0564]
 %[0.0635    0.8240    2.7811    0.0515    1.2121    7.6868]
 
%arq.x_paper= [.147, .672, 3.465, 0.566, 0.822, 8.048];
%arq.x0=[.147, .672, 3.465, 0.566, 0.822, 8.048];

 
%N/mm^2 = MPa = 1MN/(1e2cm)^2 = 1MN/1e4cm^2 = 1e2N/cm^2
 %2e5MPa = 2e7 N/cm^2

Stochast_moment = 1;
% DinamcVar = 5000;Method='MC';
% DinamcVar = 2;Method='MC';
DinamcVar = 2;Method='PC';

DinamcVar = 2;Method='R2DO';

%DinamcVar = 0;Method='FORM';
% DinamcVar = 0;Method='FERUM';


%k_p = 2000;%500
k_p = 1;
RBtype=1;%POD
%RBtype=2;%RBM
 RBtype=0;%FEM
linear_prob=1;
%RBtype=0;
 
 %% Constraint
global ResType n_sigma sign_const

%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy

arq.tpres=[0 0 0 0 0 0];% Reliability
 
%  iconst = 3;%stress
%  arq.tpres(iconst)=1;
%  ResType{iconst}='probability';
%  n_sigma(iconst)=3;
%  sign_const(iconst)=1;
%  arq.cm = -15000;
%  arq.em =  5000;
 

iconst = 1;%vol
arq.tpres(iconst)=1;
ResType{iconst}='deterministic';
n_sigma(iconst)=0;
sign_const(iconst)=1;
 arq.vol0 = 750;

%% Objective
%v_out = {vol,u,vonmis,sn,dvol,du',dsig',den'};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5 - Lambd

n_tot_obj = 5;
global mo  obj_type i_obj
mo = zeros(1,n_tot_obj*2);
obj_type=cell(1,n_tot_obj*2);

i_obj=4;% energy
arq.tobj=0;%Multi_Obj
mo([i_obj,n_tot_obj+i_obj])=1;
 obj_type{i_obj}='single';
 obj_type{n_tot_obj+i_obj}='single';
 %i_obj = 661; %specific
 arq.tobj=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random Variables

% RandVarProp (2,nrv) = [type var; which (design regions, elem, node);
% 1 type var: 1 - design variables value; 2 - material prop (E); 3 - Load
% 2 which: 1 - design regions, 2 - elem, 3 - node, 4 - Glb, 5 - factor,
% 51 - design regions

%rv_data_comp: Cell(nrv) = {components rv 1, components rv 2,...} (see which)

% RandVarDist (6,nrv) = [ distribuction; mean; std; par3; par4;....]
% 3 distribuction: 1 - Normal, 2 - Lognormal
% 4 pars: par1 - mean, par2 std.

Stochast_on=1;
nRandVar = 14;
var_tpe(1:6) = 2; % Sensibility variable E
var_tpe(7:8) = 3; % Sensibility variable Px
var_tpe(9:14) = 1; % Sensibility variable x

E=arq.E;dE=E/5;% 5,0%


%Random Variables
% 1 - E all design regions (3)
% 2 - Load (1)
% RandVarProp(1,:)=[2 2 2 3]; % type var
% RandVarProp(2,:)=[1 1 1 5]; % which type data
% rv_data_comp = {1,2,3,0};
% 
% RandVarDist(1,:)=[2 2 2 1]; % distribuction
% RandVarDist(2,:)=[E E E 1]; % Mu par1
% RandVarDist(3,:)=[dE dE dE .25]; % std par2


RandVarProp(1,:)=var_tpe; % type var
RandVarProp(2,1: 6)=51; % which type data
RandVarProp(2,7: 8)=4; % which type data 
RandVarProp(2,9:14)=5; % which type data

%(DoF node i direction x = i*3-2)
rv_data_comp = {1,2,3,4,5,6,3*3-2,6*3-2,1,2,3,4,5,6};

RandVarDist(1,1: 6)=2; % distribuction
RandVarDist(2,1: 6)=1e7/1e7; % Mu par1
RandVarDist(1,7: 8)=1; % distribuction
RandVarDist(2,7: 8)=500; % Mu par1
RandVarDist(1,9:14)=1; % distribuction
RandVarDist(2,9:14)=1; % Mu par1

RandVarDist(3,1: 5)=2e5/1e7; % std par2
RandVarDist(3,   6)=1.5e6/1e7; % std par2
RandVarDist(3,7: 8)=50; % std par2
RandVarDist(3,9:14)=0.05; % std par2

%%%%%%%%%%%%%%%%%%%%%%%%
% R. V. Elimination!! 
% R. V. Considered:
RV_consider = [5,6,11,13,14];
% RV_consider = [1:14];
nRandVar = length(RV_consider);
rv_data_comp = {rv_data_comp{RV_consider}};
RandVarDist=RandVarDist(:,RV_consider);
RandVarProp=RandVarProp(:,RV_consider);
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limit State (G=0)
G_fun = 'StructLimitState';
S_fun = 'Stoch_Analysis';
%G_const = v_out = {vol,u,vonmis,sn,Lambd,dvol,du',dsig',den',dLambd};
% 1 - vol; 2 - displ; 3 - stress; 4 - energy; 5 - Lambd

 G_const = 3;
 G_limits = 5000;
 
serieLimitState = [];

betatarg = 3.;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrv=size(RandVarProp,2);
Correl=eye(nrv);


arq=MatProp(arq);
%SIGu=arq.MatCurv(end)*.5;